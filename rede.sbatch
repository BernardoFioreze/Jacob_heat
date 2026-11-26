#!/bin/bash
#SBATCH --job-name=heat-rede
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=saidas/%x-%j.out
##SBATCH --partition=NOME_DA_PARTICAO   # descomente e ajuste se o Atlântica exigir uma partição

# Mede a velocidade de troca entre os nós e descobre qual placa de rede o Open MPI
# deve usar.
#
# Os nós do Atlântica têm 4 placas de 1 Gb/s (eno1 a eno4) em redes diferentes, e o
# Open MPI, por padrão, usa todas. O teste tem três partes:
#   1. Ping-pong entre dois processos no mesmo nó (memória compartilhada), como referência.
#   2. Ping-pong entre dois nós, com todas as placas e com cada placa sozinha:
#      latência e velocidade para mensagens de 8 B a 4 MB.
#   3. A versão MPI do solver com P=32 em 4 nós, por um número fixo de iterações,
#      com todas as placas e com cada placa sozinha: tempo por iteração.
# A placa com o menor tempo vai no paralelo.sbatch (BTL_IF=<placa>).
#
# Uso (a partir da raiz do repositório):
#   sbatch rede.sbatch [iterações do solver]      (padrão: 20000)

ITERS=${1:-20000}

set -euo pipefail
export LC_NUMERIC=C   # ponto decimal nas contas do awk, mesmo com o sistema em português
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

echo "=== Ambiente ==="
echo "Data      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job       : ${SLURM_JOB_ID:-local}"
echo "Nós       : ${SLURM_JOB_NODELIST:-$(hostname)}"
echo "MPI       : $(mpirun --version 2>/dev/null | head -1 || echo '?')"
if git rev-parse --short HEAD >/dev/null 2>&1; then
    echo "Commit    : $(git rev-parse --short HEAD)"
fi
echo
echo "Endereço de cada nó pelo nome:"
for n in $(scontrol show hostnames "${SLURM_JOB_NODELIST:-$(hostname)}"); do
    getent hosts "$n" || echo "$n: ?"
done
echo
echo "Placas de rede em $(hostname):"
ip -br addr | grep -v '^lo'
for i in /sys/class/net/*; do
    [ "$(basename "$i")" = lo ] || echo "  $(basename "$i"): $(cat "$i/speed" 2>/dev/null || echo '?') Mb/s"
done
echo

echo "=== Compilação ==="
BUILD_DIR="build/${SLURM_JOB_ID:-local-$$}"
mkdir -p "$BUILD_DIR"
trap 'rm -rf "$BUILD_DIR"' EXIT
mpicxx -O3 -march=native -std=c++17 -Wall -Wextra -o "$BUILD_DIR/heat_mpi" heat_mpi.cpp vtkWriter.cpp
mpicxx -O3 -march=native -std=c++17 -Wall -Wextra -o "$BUILD_DIR/pingpong" pingpong.cpp
echo "ok"
echo

IFACES=$(ip -br link | awk '$1 != "lo" && $2 == "UP" { print $1 }')
RESUMO=""

usar_placa() {
    if [ "$1" = todas ]; then
        unset OMPI_MCA_btl_tcp_if_include
    else
        export OMPI_MCA_btl_tcp_if_include="$1"
    fi
}

# Tira da saída do ping-pong a linha de um tamanho: "ida_us=... MBps=... Mbps=..."
medida() {
    echo "$1" | awk -v b="$2" '$1 == "PINGPONG" && $2 == "bytes=" b { print $3, $4, $5 }'
}

echo "=== 1. Ping-pong no mesmo nó (memória compartilhada) ==="
OUT=$(timeout 120 srun --nodes=1 --ntasks=2 --cpu-bind=cores "$BUILD_DIR/pingpong" 2>&1) || true
echo "$OUT" | grep -v '^PINGPONG'
RESUMO+="PINGPONG mesmo_no 8KB: $(medida "$OUT" 8192) | 4MB: $(medida "$OUT" 4194304)"$'\n'
echo

echo "=== 2. Ping-pong entre dois nós ==="
for IF in todas $IFACES; do
    echo "--- Placa: $IF"
    usar_placa "$IF"
    # Limite de 2 min: uma placa que não liga os nós pode deixar o MPI travado.
    if OUT=$(timeout 120 srun --nodes=2 --ntasks=2 --ntasks-per-node=1 --cpu-bind=cores \
                 "$BUILD_DIR/pingpong" 2>&1); then
        echo "$OUT" | grep -v '^PINGPONG'
        RESUMO+="PINGPONG placa=$IF 8KB: $(medida "$OUT" 8192) | 4MB: $(medida "$OUT" 4194304)"$'\n'
    else
        echo "$OUT" | tail -5
        RESUMO+="PINGPONG placa=$IF falhou (erro ou mais de 120 s)"$'\n'
    fi
    echo
done

echo "=== 3. Solver com P=32 em 4 nós, $ITERS iterações ==="
for IF in todas $IFACES; do
    echo "--- Placa: $IF"
    usar_placa "$IF"
    if OUT=$(timeout 180 srun --ntasks=32 --cpu-bind=cores --distribution=block:cyclic \
                 "$BUILD_DIR/heat_mpi" 1024 1e-6 "$ITERS" 2>&1); then
        echo "$OUT" | grep -E '^(Tempo laço|Cálculo|Verificação|Troca|Espera)'
        LACO=$(echo "$OUT" | awk -F' : ' '/^Tempo laço/ { print $2 }')
        US=$(awk -v l="$LACO" -v n="$ITERS" 'BEGIN { printf "%.1f", l / n * 1e6 }')
        RESUMO+="SOLVER placa=$IF laco_s=$LACO us_por_iteracao=$US"$'\n'
    else
        echo "$OUT" | tail -5
        RESUMO+="SOLVER placa=$IF falhou (erro ou mais de 180 s)"$'\n'
    fi
    echo
done

echo "=== Resumo ==="
printf '%s' "$RESUMO"
echo
echo "Referências:"
echo "  - 8 KB a 1 Gb/s levam cerca de 66 us só de transmissão; a 100 Mb/s, 655 us."
echo "  - No job 47521 (todas as placas), uma iteração do solver levou 1483.4 us."
