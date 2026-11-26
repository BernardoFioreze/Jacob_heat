#!/bin/bash
#SBATCH --job-name=heat-rede
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --exclusive
#SBATCH --time=00:20:00
#SBATCH --output=saidas/%x-%j.out
##SBATCH --partition=NOME_DA_PARTICAO   # descomente e ajuste se o Atlântica exigir uma partição

# Descobre qual placa de rede o Open MPI deve usar entre os nós.
#
# Os nós do Atlântica têm 4 placas de 1 Gb/s (eno1 a eno4) em redes diferentes, e o
# Open MPI, por padrão, usa todas. Este teste roda a versão MPI com P=32 em 4 nós por
# um número fixo de iterações: uma vez com todas as placas e uma vez com cada placa
# sozinha. A placa com o menor tempo por iteração é a que vai no paralelo.sbatch
# (BTL_IF=<placa>).
#
# Uso (a partir da raiz do repositório):
#   sbatch rede.sbatch [iterações]      (padrão: 20000, menos de 1 min por placa)

ITERS=${1:-20000}

set -euo pipefail
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
echo

echo "=== Compilação ==="
BUILD_DIR="build/${SLURM_JOB_ID:-local-$$}"
mkdir -p "$BUILD_DIR"
trap 'rm -rf "$BUILD_DIR"' EXIT
BIN="$BUILD_DIR/heat_mpi"
mpicxx -O3 -march=native -std=c++17 -Wall -Wextra -o "$BIN" heat_mpi.cpp vtkWriter.cpp
echo

RESUMO=""
IFACES=$(ip -br link | awk '$1 != "lo" && $2 == "UP" { print $1 }')

for IF in todas $IFACES; do
    echo "=== Placa: $IF ==="
    if [ "$IF" = todas ]; then
        unset OMPI_MCA_btl_tcp_if_include
    else
        export OMPI_MCA_btl_tcp_if_include="$IF"
    fi

    # Limite de 3 min: uma placa que não liga os nós pode deixar o MPI travado.
    if OUT=$(timeout 180 srun --ntasks=32 --cpu-bind=cores --distribution=block:cyclic \
                 "$BIN" 1024 1e-6 "$ITERS" 2>&1); then
        echo "$OUT" | grep -E '^(Iterações|Tempo laço|Cálculo|Verificação|Troca|Espera)'
        LACO=$(echo "$OUT" | awk -F' : ' '/^Tempo laço/ { print $2 }')
        US=$(awk -v l="$LACO" -v n="$ITERS" 'BEGIN { printf "%.1f", l / n * 1e6 }')
        LINHA="REDE placa=$IF iteracoes=$ITERS laco_s=$LACO us_por_iteracao=$US"
        echo "$LINHA"
        RESUMO+="$LINHA"$'\n'
    else
        echo "$OUT" | tail -5
        LINHA="REDE placa=$IF falhou (erro ou mais de 180 s)"
        echo "$LINHA"
        RESUMO+="$LINHA"$'\n'
    fi
    echo
done

echo "=== Resumo (menor us_por_iteracao = melhor placa) ==="
printf '%s' "$RESUMO"
echo "Referência: no job 47521 (todas as placas), uma iteração levou 1483.4 us."
