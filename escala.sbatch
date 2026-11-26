#!/bin/bash
#SBATCH --job-name=heat-escala
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --ntasks-per-core=2
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=saidas/%x-%j.out
##SBATCH --partition=NOME_DA_PARTICAO   # descomente e ajuste se o Atlântica exigir uma partição

# Prévia rápida da escalabilidade, antes dos jobs longos.
#
# Num único job de 4 nós, roda a versão MPI com P = 2, 4, 8, 16, 32 e 64 por um número
# fixo de iterações, cada P com os nós do enunciado (1 nó até 8 processos, 2 com 16,
# 4 com 32 e 64, este com HT). Para cada P, mostra:
#   - o tempo de cada fase, por processo (--por-processo), com o nó, o soquete e o núcleo;
#   - o tempo por iteração e a estimativa de uma repetição completa (2 449 575 iterações),
#     do speed-up e da eficiência, comparando com o T1 medido (16 515,47 s, job 47421).
#
# Uso (a partir da raiz do repositório):
#   sbatch escala.sbatch [iterações]      (padrão: 20000, uns 5 min no total)
# Variável opcional: BTL_IF (placa de rede entre os nós, padrão eno1; vazio = todas).

ITERS=${1:-20000}
T1=16515.47          # tempo do laço do T1 (job 47421), em s
FULL_ITERS=2449575   # iterações até convergir com N=1024, tol=1e-6

set -euo pipefail
export LC_NUMERIC=C   # ponto decimal nas contas do awk, mesmo com o sistema em português
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

BTL_IF=${BTL_IF-eno1}
if [ -n "$BTL_IF" ]; then
    export OMPI_MCA_btl_tcp_if_include="$BTL_IF"
fi

echo "=== Ambiente ==="
echo "Data      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job       : ${SLURM_JOB_ID:-local}"
echo "Nós       : ${SLURM_JOB_NODELIST:-$(hostname)}"
echo "CPU       : $(lscpu | sed -n 's/^Model name:[[:space:]]*//p;s/^Nome do modelo:[[:space:]]*//p' | head -1)"
echo "MPI       : $(mpirun --version 2>/dev/null | head -1 || echo '?')"
echo "Rede MPI  : ${BTL_IF:-todas as placas}"
if git rev-parse --short HEAD >/dev/null 2>&1; then
    echo "Commit    : $(git rev-parse --short HEAD)"
fi
echo "Iterações : $ITERS por P"
echo

echo "=== Compilação ==="
BUILD_DIR="build/${SLURM_JOB_ID:-local-$$}"
mkdir -p "$BUILD_DIR"
trap 'rm -rf "$BUILD_DIR"' EXIT
BIN="$BUILD_DIR/heat_mpi"
mpicxx -O3 -march=native -std=c++17 -Wall -Wextra -o "$BIN" heat_mpi.cpp vtkWriter.cpp
echo "ok"
echo

RESUMO=""

# P, nós, processos por nó, ligação (cores = 1 por núcleo físico; threads = HT)
for CFG in "2 1 2 cores" "4 1 4 cores" "8 1 8 cores" "16 2 8 cores" "32 4 8 cores" "64 4 16 threads"; do
    read -r P NODES PER_NODE BIND <<< "$CFG"
    if [ "$BIND" = cores ]; then PER_CORE=1; else PER_CORE=2; fi
    echo "=== P=$P ($NODES nó(s), $PER_NODE por nó, ligação por $BIND) ==="

    if OUT=$(timeout 300 srun --nodes="$NODES" --ntasks="$P" --ntasks-per-node="$PER_NODE" \
                 --ntasks-per-core="$PER_CORE" --cpu-bind="$BIND" --distribution=block:cyclic \
                 "$BIN" 1024 1e-6 "$ITERS" --por-processo 2>&1); then
        echo "$OUT" | grep -E '^(Tempo laço|Fases|Cálculo|Verificação|Troca|Espera|Por processo|PROC)'
        LACO=$(echo "$OUT" | awk -F' : ' '/^Tempo laço/ { print $2 }')
        LINHA=$(awk -v p="$P" -v l="$LACO" -v n="$ITERS" -v full="$FULL_ITERS" -v t1="$T1" 'BEGIN {
            us  = l / n * 1e6
            rep = us * 1e-6 * full
            s   = t1 / rep
            printf "ESCALA P=%d us_por_iteracao=%.1f repeticao_min=%.1f speedup=%.2f eficiencia=%.2f", p, us, rep / 60, s, s / p
        }')
        echo "$LINHA"
        RESUMO+="$LINHA"$'\n'
    else
        echo "$OUT" | tail -8
        RESUMO+="ESCALA P=$P falhou (erro ou mais de 300 s)"$'\n'
    fi
    echo
done

echo "=== Resumo (estimativas para N=1024, tol=1e-6, a partir de $ITERS iterações) ==="
printf '%s' "$RESUMO"
