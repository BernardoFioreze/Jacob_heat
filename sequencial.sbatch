#!/bin/bash
#SBATCH --job-name=heat-seq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --time=08:00:00
#SBATCH --output=saidas/%x-%j.out
##SBATCH --partition=NOME_DA_PARTICAO   # descomente e ajuste se o Atlântica exigir uma partição

# Versão sequencial do solver de calor 2D (Jacobi) no Atlântica.
#
# Uso (a partir da raiz do repositório):
#   sbatch sequencial.sbatch [N] [tol] [repetições]
# Padrão: N=1024, tol=1e-6, 1 repetição (cerca de 80 min num i7-1355U, mas
# mais de 3 h no Xeon E5520 do Atlântica: o job 47341 foi cancelado nas 3 h).
# Limite de 8 h; o cluster permite até 16 h para jobs de 1 nó, mas só 2 jobs
# na fila por usuário.
#
# A saída vai para saidas/heat-seq-<jobid>.out.

N=${1:-1024}
TOL=${2:-1e-6}
REPS=${3:-1}

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

# Se o compilador vier de um módulo no Atlântica, carregue aqui, por exemplo:
# module load gcc

echo "=== Ambiente ==="
echo "Data      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job       : ${SLURM_JOB_ID:-local}"
echo "Nó        : $(hostname)"
echo "CPU       : $(lscpu | sed -n 's/^Model name:[[:space:]]*//p;s/^Nome do modelo:[[:space:]]*//p' | head -1)"
echo "Compilador: $(${CXX:-g++} --version | head -1)"
if git rev-parse --short HEAD >/dev/null 2>&1; then
    echo "Commit    : $(git rev-parse --short HEAD)"
fi
echo "Parâmetros: N=$N tol=$TOL repetições=$REPS"
echo

# Compila no próprio nó de cálculo, porque -march=native gera código para a CPU de onde compila.
# O binário vai para uma pasta temporária do job, para jobs simultâneos não se atrapalharem.
echo "=== Compilação ==="
BUILD_DIR=$(mktemp -d)
trap 'rm -rf "$BUILD_DIR"' EXIT
BIN="$BUILD_DIR/heat"
make --always-make BIN="$BIN"
echo

# Dentro do Slurm, srun prende o processo a um núcleo; fora dele, roda direto.
RUN=""
if [ -n "${SLURM_JOB_ID:-}" ] && command -v srun >/dev/null 2>&1; then
    RUN="srun --ntasks=1 --cpu-bind=cores"
fi

for rep in $(seq 1 "$REPS"); do
    echo "=== Repetição $rep de $REPS ==="
    OUT=$($RUN "$BIN" "$N" "$TOL")
    echo "$OUT"

    # Uma linha por repetição, fácil de filtrar com: grep RESULTADO saidas/*.out
    echo "$OUT" | awk -v rep="$rep" -F' : |: ' '
        /^N /                  { n = $2 }
        /^Tolerância/          { tol = $2 }
        /^Iterações/           { it = $2 }
        /^Soma de T/           { soma = $2 }
        /^Convergiu/           { conv = $2 }
        /^Tempo laço/          { laco = $2 }
        /^Tempo total/         { total = $2 }
        /^Fração fora do laço/ { frac = $2 }
        END { printf "RESULTADO rep=%s N=%s tol=%s iteracoes=%s soma=%s convergiu=%s laco_s=%s total_s=%s fracao_seq=%s\n",
                     rep, n, tol, it, soma, conv, laco, total, frac }'
    echo
done
