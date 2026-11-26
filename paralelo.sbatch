#!/bin/bash
#SBATCH --job-name=heat-mpi
#SBATCH --exclusive
#SBATCH --output=saidas/%x-%j.out
##SBATCH --partition=NOME_DA_PARTICAO   # descomente e ajuste se o Atlântica exigir uma partição

# Versão paralela (MPI) do solver de calor 2D no Atlântica.
#
# Uso (a partir da raiz do repositório):
#   sbatch -J <nome> --nodes=<nós> --ntasks-per-node=<P/nós> --time=<limite> \
#          paralelo.sbatch P [repetições] [primeira] [N] [tol]
# Padrão: 10 repetições, começando na 1, N=1024, tol=1e-6.
#
# Nós, tarefas e tempo vão na linha de comando, porque mudam com P:
#   sbatch -J heat-p2  --nodes=1 --ntasks-per-node=2  --time=16:00:00 paralelo.sbatch 2 5 1
#   sbatch -J heat-p2  --nodes=1 --ntasks-per-node=2  --time=16:00:00 paralelo.sbatch 2 5 6
#   sbatch -J heat-p4  --nodes=1 --ntasks-per-node=4  --time=16:00:00 paralelo.sbatch 4 5 1
#   sbatch -J heat-p4  --nodes=1 --ntasks-per-node=4  --time=16:00:00 paralelo.sbatch 4 5 6
#   sbatch -J heat-p8  --nodes=1 --ntasks-per-node=8  --time=16:00:00 paralelo.sbatch 8
#   sbatch -J heat-p16 --nodes=2 --ntasks-per-node=8  --time=08:00:00 paralelo.sbatch 16
#   sbatch -J heat-p32 --nodes=4 --ntasks-per-node=8  --time=04:00:00 paralelo.sbatch 32
#   sbatch -J heat-p64 --nodes=4 --ntasks-per-node=16 --ntasks-per-core=2 --time=04:00:00 paralelo.sbatch 64
#
# "primeira" é o número da primeira repetição do job. Serve para dividir as 10
# repetições em mais de um job (1 a 5 num, 6 a 10 no outro) ou para completar
# as que faltaram se um job for cancelado por tempo.
#
# Variáveis de ambiente opcionais:
#   MPI_SRC   fontes da versão MPI (padrão: "heat_mpi.cpp vtkWriter.cpp")
#   MPICXX    compilador MPI (padrão: mpicxx)
#   LAUNCHER  srun (padrão), mpirun, ou um comando completo
#   BTL_IF    placa de rede que o Open MPI usa entre os nós (ex.: BTL_IF=eno1);
#             vazio usa todas. Descubra a melhor com: sbatch rede.sbatch
#
# A saída vai para saidas/<nome>-<jobid>.out.

P=${1:?informe P (número de processos)}
REPS=${2:-10}
FIRST=${3:-1}
N=${4:-1024}
TOL=${5:-1e-6}
MPI_SRC=${MPI_SRC:-heat_mpi.cpp vtkWriter.cpp}
MPICXX=${MPICXX:-mpicxx}
CXXFLAGS=${CXXFLAGS:--O3 -march=native -std=c++17 -Wall -Wextra}
LAUNCHER=${LAUNCHER:-srun}
BTL_IF=${BTL_IF:-}
if [ -n "$BTL_IF" ]; then
    export OMPI_MCA_btl_tcp_if_include="$BTL_IF"
fi

# Iterações e soma de T do T1 no Atlântica (job 47421), para conferir a equivalência.
REF_ITER=""
REF_SOMA=""
if [ "$N" = 1024 ] && [ "$TOL" = 1e-6 ]; then
    REF_ITER=2449575
    REF_SOMA=167156330.0959675014
fi

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"

# Se o MPI vier de um módulo no Atlântica, carregue aqui, por exemplo:
# module load openmpi

# Núcleos físicos do nó (sem contar HT) e processos por nó.
CORES=$(lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l)
NNODES=${SLURM_NNODES:-1}
PER_NODE=$(( (P + NNODES - 1) / NNODES ))

# Um processo por núcleo físico; só com mais processos que núcleos (P=64) usa HT.
if [ "$PER_NODE" -gt "$CORES" ]; then BIND=threads; else BIND=cores; fi

echo "=== Ambiente ==="
echo "Data      : $(date '+%Y-%m-%d %H:%M:%S')"
echo "Job       : ${SLURM_JOB_ID:-local}"
echo "Nós       : ${SLURM_JOB_NODELIST:-$(hostname)} ($NNODES)"
echo "CPU       : $(lscpu | sed -n 's/^Model name:[[:space:]]*//p;s/^Nome do modelo:[[:space:]]*//p' | head -1)"
echo "Núcleos   : $CORES físicos por nó, $PER_NODE processos por nó, ligação por $BIND"
echo "Compilador: $($MPICXX --version | head -1)"
echo "MPI       : $(mpirun --version 2>/dev/null | head -1 || echo '?')"
echo "Rede MPI  : ${BTL_IF:-todas as placas}"
if git rev-parse --short HEAD >/dev/null 2>&1; then
    echo "Commit    : $(git rev-parse --short HEAD)"
fi
echo "Parâmetros: P=$P N=$N tol=$TOL repetições=$FIRST a $((FIRST + REPS - 1))"
echo

# Confere se o job pediu os nós e tarefas certos para este P.
if [ -n "${SLURM_JOB_ID:-}" ]; then
    TASKS=${SLURM_NTASKS:-}
    if [ -z "$TASKS" ]; then
        PN=${SLURM_NTASKS_PER_NODE:-0}   # pode vir como "8(x2)"
        PN=${PN%%[(]*}
        TASKS=$(( NNODES * PN ))
    fi
    if [ "$TASKS" -ne "$P" ]; then
        echo "ERRO: o job tem $TASKS tarefas, mas P=$P. Confira --nodes e --ntasks-per-node."
        exit 1
    fi
    if   [ "$P" -le 8 ];  then WANT=1
    elif [ "$P" -le 16 ]; then WANT=2
    else                       WANT=4
    fi
    if [ "$NNODES" -ne "$WANT" ]; then
        echo "AVISO: o enunciado pede $WANT nó(s) para P=$P, e este job tem $NNODES."
        echo
    fi
fi

# Compila no nó de cálculo (-march=native). O binário fica numa pasta do job dentro
# do repositório, porque precisa estar visível em todos os nós; /tmp é local de cada nó.
echo "=== Compilação ==="
read -ra SOURCES <<< "$MPI_SRC"
for src in "${SOURCES[@]}"; do
    if [ ! -f "$src" ]; then
        echo "ERRO: fonte $src não encontrada. Ajuste MPI_SRC (ex.: MPI_SRC=meu_arquivo.cpp sbatch ...)."
        exit 1
    fi
done
BUILD_DIR="build/${SLURM_JOB_ID:-local-$$}"
mkdir -p "$BUILD_DIR"
trap 'rm -rf "$BUILD_DIR"' EXIT
BIN="$BUILD_DIR/heat_mpi"
# shellcheck disable=SC2086
echo "$MPICXX $CXXFLAGS -o $BIN ${SOURCES[*]}"
# shellcheck disable=SC2086
$MPICXX $CXXFLAGS -o "$BIN" "${SOURCES[@]}"
echo

# Como lançar os P processos.
case "$LAUNCHER" in
    srun)
        # block:cyclic alterna os processos entre os dois soquetes do nó,
        # para usar os dois controladores de memória.
        RUN=(srun --ntasks="$P" --cpu-bind="$BIND" --distribution=block:cyclic)
        ;;
    mpirun)
        if [ "$BIND" = threads ]; then
            RUN=(mpirun -np "$P" --use-hwthread-cpus --map-by socket --bind-to hwthread)
        else
            RUN=(mpirun -np "$P" --map-by socket --bind-to core)
        fi
        ;;
    *)
        read -ra RUN <<< "$LAUNCHER"
        ;;
esac
echo "Lançador  : ${RUN[*]:-(direto)}"
echo

LAST=$((FIRST + REPS - 1))
for rep in $(seq "$FIRST" "$LAST"); do
    echo "=== Repetição $rep (P=$P) ==="
    if ! OUT=$("${RUN[@]}" "$BIN" "$N" "$TOL"); then
        echo "$OUT"
        echo "ERRO: a repetição $rep falhou."
        echo
        continue
    fi
    echo "$OUT"

    # Uma linha por repetição, fácil de filtrar com: grep RESULTADO saidas/*.out
    # equivalente: mesmo número de iterações do T1; soma_igual: mesma soma de T, dígito por dígito.
    # Fases: média entre os processos (calc, verif, troca) e espera estimada pelo mais lento.
    echo "$OUT" | awk -v rep="$rep" -v p="$P" -v ref="$REF_ITER" -v refsoma="$REF_SOMA" -F' : |: ' '
        /^N /                  { n = $2 }
        /^Tolerância/          { tol = $2 }
        /^Iterações/           { it = $2 }
        /^Soma de T/           { soma = $2 }
        /^Convergiu/           { conv = $2 }
        /^Tempo laço/          { laco = $2 }
        /^Tempo total/         { total = $2 }
        /^Fração fora do laço/ { frac = $2 }
        /^Cálculo/             { split($2, v, " "); calc = v[2] }
        /^Verificação/         { split($2, v, " "); verif = v[2] }
        /^Troca/               { split($2, v, " "); troca = v[2] }
        /^Espera verif/        { espv = $2 }
        /^Espera troca/        { espt = $2 }
        END {
            eq = (ref == "")     ? "?" : ((it == ref)       ? "sim" : "NAO")
            sq = (refsoma == "") ? "?" : (((soma "") == (refsoma "")) ? "sim" : "NAO")
            printf "RESULTADO P=%s rep=%s N=%s tol=%s iteracoes=%s soma=%s convergiu=%s equivalente=%s soma_igual=%s laco_s=%s total_s=%s fracao_seq=%s calc_s=%s verif_s=%s troca_s=%s espera_verif_s=%s espera_troca_s=%s\n",
                   p, rep, n, tol, it, soma, conv, eq, sq, laco, total, frac, calc, verif, troca, espv, espt
        }'
    echo
done
