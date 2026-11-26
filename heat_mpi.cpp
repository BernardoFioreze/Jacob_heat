// Condução de calor 2D em regime permanente numa placa de 4 materiais.
// Diferenças finitas + método de Jacobi, versão paralela com MPI (modelo de fases paralelas).
//
// Cada processo fica com um bloco de linhas consecutivas da grade, mais uma linha
// fantasma em cima e outra embaixo. A cada iteração, todos repetem três fases:
//   1. Cálculo: varredura de Jacobi no próprio bloco, com as mesmas contas da versão sequencial.
//   2. Verificação: MPI_Allreduce da soma de |novo - antigo|; todos recebem o mesmo total
//      e decidem juntos se param.
//   3. Troca: cada processo manda a primeira e a última linha do bloco aos vizinhos de cima
//      e de baixo (MPI_Isend/MPI_Irecv nas duas direções ao mesmo tempo, depois MPI_Waitall),
//      que as guardam nas linhas fantasma para a próxima iteração.
//
// Cada célula é calculada com os mesmos valores e na mesma ordem de operações da versão
// sequencial, então o campo final é idêntico bit a bit, e o número de iterações também.
//
// Uso: mpirun -np P ./heat_mpi [N] [tol] [maxIter] [--vtk] [--print]
//   (mesmos argumentos da versão sequencial; precisa de N >= P)

// Só a interface em C do MPI; a de C++ (obsoleta) gera avisos de compilação no Open MPI.
#define OMPI_SKIP_MPICXX 1
#define MPICH_SKIP_MPICXX 1
#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>
#include "solverParameters.h"
#include "vtkWriter.h"

struct Settings {
    int    N       = default_N;
    double tol     = default_tol;
    long   maxIter = default_maxIter;
    bool   vtk     = false;
    bool   print   = false;
};

// Tempo acumulado em cada fase, em um processo.
struct PhaseTimes {
    double calc  = 0.0;  // fase 1: condições de contorno + varredura
    double verif = 0.0;  // fase 2: MPI_Allreduce do critério de parada
    double troca = 0.0;  // fase 3: MPI_Sendrecv das linhas de borda
};

// Material de uma célula pela posição global (i, j), contando as células fantasma.
// Mesma regra da versão sequencial: quatro quadrantes.
static int materialAt(int i, int j, int N) {
    if      (i <= N / 2 && j <= N / 2) return 1;
    else if (i <= N / 2)               return 2;
    else if (j <= N / 2)               return 3;
    else                               return 4;
}

// Condutividade na face entre duas células: média harmônica se os materiais forem diferentes.
static double faceConductivity(int mCenter, int mNeighbor) {
    double kc = conductivity[mCenter];
    double kn = conductivity[mNeighbor];
    return (mCenter != mNeighbor) ? 2.0 * kc * kn / (kc + kn) : kc;
}

// Divide as N linhas do interior entre P processos; os primeiros N % P ficam com uma a mais.
// first é a linha global (1..N) onde o bloco começa; count é quantas linhas ele tem.
static void rowBlock(int N, int P, int rank, int& first, int& count) {
    int base = N / P, extra = N % P;
    count = base + (rank < extra ? 1 : 0);
    first = 1 + rank * base + (rank < extra ? rank : extra);
}

class ParallelHeatSolver {
public:
    ParallelHeatSolver(const Settings& s, int rank, int size)
        : N(s.N), W(s.N + 2),
          dx(xLength / (s.N + 1)), dy(yLength / (s.N + 1)),
          tol(s.tol), maxIter(s.maxIter),
          rank(rank), size(size),
          up(rank > 0 ? rank - 1 : MPI_PROC_NULL),
          down(rank < size - 1 ? rank + 1 : MPI_PROC_NULL) {
        rowBlock(N, size, rank, first, L);
        long cells = static_cast<long>(L + 2) * W;
        cE.assign(cells, 0.0); cW.assign(cells, 0.0);
        cN.assign(cells, 0.0); cS.assign(cells, 0.0);
        robinA.assign(L + 2, 0.0); robinB.assign(L + 2, 0.0);
        bufferA.assign(cells, T_initial);
        bufferB.assign(cells, T_initial);
        T    = bufferA.data();
        Told = bufferB.data();
        buildCoefficients();
        setDirichletRows(T);
        setDirichletRows(Told);
    }

    // Laço de Jacobi nas três fases. Retorna o número de iterações feitas.
    long solve(PhaseTimes& pt) {
        long iteration = 0;
        avgChange = 1.0;
        const double totalCells = static_cast<double>(N) * N;

        while (avgChange > tol && iteration < maxIter) {
            // O campo da iteração anterior passa a ser o "antigo" (troca de ponteiros, sem cópia).
            // As linhas fantasma dele já chegaram na fase 3 da iteração anterior.
            std::swap(T, Told);

            // Fase 1: cálculo local.
            double t0 = MPI_Wtime();
            applyBoundaryConditions(Told);
            double localChange = jacobiSweep();
            double t1 = MPI_Wtime();

            // Fase 2: verificação da parada. Todos recebem a mesma soma global,
            // então todos tomam a mesma decisão na mesma iteração.
            double globalChange = 0.0;
            MPI_Allreduce(&localChange, &globalChange, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            avgChange = globalChange / totalCells;
            ++iteration;
            double t2 = MPI_Wtime();

            pt.calc  += t1 - t0;
            pt.verif += t2 - t1;

            // Fase 3: troca das bordas, só se houver uma próxima iteração.
            if (avgChange > tol && iteration < maxIter) {
                exchangeHalos(T);
                pt.troca += MPI_Wtime() - t2;
            }
        }
        return iteration;
    }

    double residual() const { return avgChange; }

    // Junta o campo inteiro no processo 0 (linhas 1..N, com as colunas fantasma).
    // Nos outros processos, field fica vazio.
    void gatherField(std::vector<double>& field) const {
        std::vector<int> counts, displs;
        if (rank == 0) {
            field.assign(static_cast<long>(W) * W, 0.0);
            counts.resize(size);
            displs.resize(size);
            for (int r = 0; r < size; ++r) {
                int f, c;
                rowBlock(N, size, r, f, c);
                counts[r] = c * W;
                displs[r] = f * W;
            }
        }
        MPI_Gatherv(T + idx(1, 0), L * W, MPI_DOUBLE,
                    rank == 0 ? field.data() : nullptr,
                    counts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

private:
    const int N;       // células do interior por direção (grade inteira)
    const int W;       // posições por linha no array (N + 2, com as células fantasma)
    const double dx, dy;
    const double tol;
    const long maxIter;
    const int rank, size;
    const int up, down;  // vizinhos de cima e de baixo (MPI_PROC_NULL nas pontas)
    int first = 0;       // linha global da primeira linha do bloco
    int L = 0;           // linhas do bloco (sem as fantasma)

    // Arrays locais de (L+2) x W posições: linha 0 e linha L+1 são as fantasma.
    std::vector<double> cE, cW, cN, cS;   // coeficientes já divididos por a_P
    std::vector<double> robinA, robinB;   // borda direita: T[i][N+1] = robinA[i]*T[i][N] + robinB[i]
    std::vector<double> bufferA, bufferB; // os dois campos de temperatura
    double* T    = nullptr;  // campo novo (escrito nesta iteração)
    double* Told = nullptr;  // campo antigo (lido nesta iteração)
    double avgChange = 0.0;

    long idx(int i, int j) const { return static_cast<long>(i) * W + j; }
    int globalRow(int li) const { return first + li - 1; }

    // Mesmas contas da versão sequencial, para a linha global de cada linha local.
    void buildCoefficients() {
        for (int li = 1; li <= L; ++li) {
            int i = globalRow(li);
            for (int j = 1; j <= N; ++j) {
                int m = materialAt(i, j, N);
                double aE = faceConductivity(m, materialAt(i, j + 1, N)) * dy / dx;
                double aW = faceConductivity(m, materialAt(i, j - 1, N)) * dy / dx;
                double aN = faceConductivity(m, materialAt(i - 1, j, N)) * dx / dy;
                double aS = faceConductivity(m, materialAt(i + 1, j, N)) * dx / dy;
                double aP = aE + aW + aN + aS;
                cE[idx(li, j)] = aE / aP;
                cW[idx(li, j)] = aW / aP;
                cN[idx(li, j)] = aN / aP;
                cS[idx(li, j)] = aS / aP;
            }
            double k = conductivity[materialAt(i, N, N)];
            robinA[li] = k / (k + h * dx);
            robinB[li] = h * dx * Tf / (k + h * dx);
        }
    }

    // Só o primeiro processo tem a borda superior, e só o último tem a inferior.
    void setDirichletRows(double* t) {
        for (int j = 0; j < W; ++j) {
            if (rank == 0)        t[idx(0, j)]     = T_top;
            if (rank == size - 1) t[idx(L + 1, j)] = T_bottom;
        }
    }

    void applyBoundaryConditions(double* t) {
        for (int li = 1; li <= L; ++li) {
            double* row = t + idx(li, 0);
            row[0]     = row[1];                             // esquerda: Neumann (adiabática)
            row[N + 1] = robinA[li] * row[N] + robinB[li];   // direita: Robin (convecção)
        }
    }

    // Igual à da versão sequencial, só que nas L linhas do bloco.
    double jacobiSweep() {
        double changeSum = 0.0;
        for (int li = 1; li <= L; ++li) {
            const double* upRow   = Told + idx(li - 1, 0);
            const double* center  = Told + idx(li, 0);
            const double* downRow = Told + idx(li + 1, 0);
            const double* e = &cE[idx(li, 0)];
            const double* w = &cW[idx(li, 0)];
            const double* n = &cN[idx(li, 0)];
            const double* s = &cS[idx(li, 0)];
            double* out = T + idx(li, 0);

            for (int j = 1; j <= N; ++j) {
                double v = e[j] * center[j + 1] + w[j] * center[j - 1]
                         + n[j] * upRow[j]      + s[j] * downRow[j];
                changeSum += std::fabs(v - center[j]);
                out[j] = v;
            }
        }
        return changeSum;
    }

    // A primeira linha do bloco vai para o vizinho de cima, que a guarda na sua linha
    // fantasma de baixo; a última vai para o de baixo, que a guarda na fantasma de cima.
    // As quatro operações (duas recepções e dois envios) começam juntas e só então o
    // processo espera todas terminarem: as mensagens das duas direções viajam ao mesmo
    // tempo, e a troca custa cerca de uma latência de rede, não duas.
    // Nas pontas, MPI_PROC_NULL faz a operação não enviar nem receber nada, e as linhas
    // Dirichlet ficam como estão.
    void exchangeHalos(double* t) {
        MPI_Request req[4];
        MPI_Irecv(t + idx(0, 1),     N, MPI_DOUBLE, up,   1, MPI_COMM_WORLD, &req[0]);  // última linha do de cima
        MPI_Irecv(t + idx(L + 1, 1), N, MPI_DOUBLE, down, 0, MPI_COMM_WORLD, &req[1]);  // primeira do de baixo
        MPI_Isend(t + idx(1, 1),     N, MPI_DOUBLE, up,   0, MPI_COMM_WORLD, &req[2]);  // minha primeira, para cima
        MPI_Isend(t + idx(L, 1),     N, MPI_DOUBLE, down, 1, MPI_COMM_WORLD, &req[3]);  // minha última, para baixo
        MPI_Waitall(4, req, MPI_STATUSES_IGNORE);
    }
};

// Soma das temperaturas do interior, na mesma ordem da versão sequencial
// (por isso dá o mesmo valor, até o último dígito).
static double temperatureSum(const std::vector<double>& field, int N) {
    int W = N + 2;
    double sum = 0.0;
    for (int i = 1; i <= N; ++i) {
        const double* row = field.data() + static_cast<long>(i) * W;
        for (int j = 1; j <= N; ++j) sum += row[j];
    }
    return sum;
}

static void printField(const std::vector<double>& field, int N) {
    int W = N + 2;
    std::printf("-----------Solução-----------\n");
    for (int i = 1; i <= N; ++i) {
        for (int j = 1; j <= N; ++j) std::printf("%.2f\t", field[static_cast<long>(i) * W + j]);
        std::printf("\n");
    }
}

static void writeVTK(const std::vector<double>& field, int N) {
    int W = N + 2;
    double dx = xLength / (N + 1), dy = yLength / (N + 1);
    exportToVTK(field.data(), N, W, dx, dy, "temperatureField.vtk");

    std::vector<double> kField(static_cast<long>(W) * W);
    for (int i = 0; i < W; ++i)
        for (int j = 0; j < W; ++j)
            kField[static_cast<long>(i) * W + j] = conductivity[materialAt(i, j, N)];
    exportToVTK(kField.data(), N, W, dx, dy, "K_conductivity.vtk");
}

// Encerra todos os processos com uma mensagem (impressa só pelo processo 0).
static void fail(int rank, const char* msg, const char* arg = "") {
    if (rank == 0) std::fprintf(stderr, msg, arg);
    MPI_Finalize();
    std::exit(1);
}

static Settings parseArgs(int argc, char** argv, int rank, int size) {
    Settings s;
    int positional = 0;
    for (int a = 1; a < argc; ++a) {
        if      (std::strcmp(argv[a], "--vtk") == 0)   s.vtk = true;
        else if (std::strcmp(argv[a], "--print") == 0) s.print = true;
        else if (positional == 0) { s.N = std::atoi(argv[a]);        ++positional; }
        else if (positional == 1) { s.tol = std::atof(argv[a]);      ++positional; }
        else if (positional == 2) { s.maxIter = std::atol(argv[a]);  ++positional; }
        else fail(rank, "Argumento desconhecido: %s\n", argv[a]);
    }
    if (s.N < 2 || s.tol <= 0.0 || s.maxIter < 1)
        fail(rank, "Uso: %s [N>=2] [tol>0] [maxIter>=1] [--vtk] [--print]\n", argv[0]);
    if (s.N < size)
        fail(rank, "N precisa ser pelo menos o número de processos (cada um fica com uma linha ou mais).\n");
    return s;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Settings settings = parseArgs(argc, argv, rank, size);

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    ParallelHeatSolver solver(settings, rank, size);   // alocação e coeficientes do bloco
    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    PhaseTimes pt;
    long iterations = solver.solve(pt);                 // laço iterativo (três fases)
    MPI_Barrier(MPI_COMM_WORLD);
    double t2 = MPI_Wtime();

    std::vector<double> field;
    solver.gatherField(field);                          // campo inteiro no processo 0
    double sum = 0.0;
    if (rank == 0) {
        sum = temperatureSum(field, settings.N);
        if (settings.print) printField(field, settings.N);
        if (settings.vtk)   writeVTK(field, settings.N);
    }
    double t3 = MPI_Wtime();

    // Tempo de cada fase: mínimo, média e máximo entre os processos.
    double local[3] = {pt.calc, pt.verif, pt.troca};
    double tmin[3], tsum[3], tmax[3];
    MPI_Reduce(local, tmin, 3, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(local, tsum, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(local, tmax, 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double tInit = t1 - t0, tSolve = t2 - t1, tOut = t3 - t2, tTotal = t3 - t0;
        double tavg[3];
        for (int k = 0; k < 3; ++k) tavg[k] = tsum[k] / size;

        std::printf("N               : %d\n", settings.N);
        std::printf("Processos       : %d\n", size);
        std::printf("Tolerância      : %g\n", settings.tol);
        std::printf("Iterações       : %ld\n", iterations);
        std::printf("Resíduo final   : %.10e\n", solver.residual());
        std::printf("Soma de T       : %.10f\n", sum);
        std::printf("Convergiu       : %s\n", solver.residual() <= settings.tol ? "sim" : "NÃO (atingiu maxIter)");
        std::printf("Tempo init  (s) : %.6f\n", tInit);
        std::printf("Tempo laço  (s) : %.6f\n", tSolve);
        std::printf("Tempo saída (s) : %.6f\n", tOut);
        std::printf("Tempo total (s) : %.6f\n", tTotal);
        std::printf("Fração fora do laço: %.6f\n", (tInit + tOut) / tTotal);
        std::printf("Fases (s), entre os processos: mín méd máx\n");
        std::printf("Cálculo     (s) : %.6f %.6f %.6f\n", tmin[0], tavg[0], tmax[0]);
        std::printf("Verificação (s) : %.6f %.6f %.6f\n", tmin[1], tavg[1], tmax[1]);
        std::printf("Troca       (s) : %.6f %.6f %.6f\n", tmin[2], tavg[2], tmax[2]);
        // O processo que chega por último a uma comunicação quase não espera, então o
        // mínimo entre os processos estima o custo da comunicação em si, e média - mínimo
        // estima quanto, em média, cada processo ficou esperando o mais lento.
        std::printf("Espera verif.   (s) : %.6f\n", tavg[1] - tmin[1]);
        std::printf("Espera troca    (s) : %.6f\n", tavg[2] - tmin[2]);
    }

    MPI_Finalize();
    return 0;
}
