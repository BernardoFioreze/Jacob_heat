// Condução de calor 2D em regime permanente numa placa de 4 materiais.
// Diferenças finitas + método de Jacobi, versão sequencial.
//
// Baseado em: https://github.com/thermeng/2D-Heat-Conduction-Solver-for-CompositePlate
//
// Uso: ./heat [N] [tol] [maxIter] [--vtk] [--print]
//   N       : células do interior em cada direção (padrão 16)
//   tol     : para quando a variação média por célula fica <= tol (padrão 1e-6)
//   maxIter : trava de segurança no número de iterações (padrão 1e8)
//   --vtk   : escreve temperatureField.vtk e K_conductivity.vtk
//   --print : imprime o campo de temperatura no terminal (só para N pequeno)

#include <chrono>
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

class HeatConductionSolver {
public:
    explicit HeatConductionSolver(const Settings& s)
        : N(s.N), W(s.N + 2),
          dx(xLength / (s.N + 1)), dy(yLength / (s.N + 1)),
          tol(s.tol), maxIter(s.maxIter),
          material(cells()), cE(cells()), cW(cells()), cN(cells()), cS(cells()),
          robinA(W), robinB(W),
          bufferA(cells(), T_initial), bufferB(cells(), T_initial) {
        T    = bufferA.data();
        Told = bufferB.data();
        buildMaterialMap();
        buildCoefficients();
        setDirichletRows(T);
        setDirichletRows(Told);
    }

    // Laço de Jacobi. Retorna o número de iterações feitas.
    long solve() {
        long iteration = 0;
        avgChange = 1.0;

        while (avgChange > tol && iteration < maxIter) {
            // O campo da iteração anterior passa a ser o "antigo" (troca de ponteiros, sem cópia).
            std::swap(T, Told);

            // Células fantasma das bordas esquerda e direita, calculadas a partir do campo antigo.
            applyBoundaryConditions(Told);

            // Atualiza o interior e devolve a soma de |novo - antigo|.
            double changeSum = jacobiSweep();

            avgChange = changeSum / (static_cast<double>(N) * N);
            ++iteration;
        }
        return iteration;
    }

    double residual() const { return avgChange; }

    // Soma das temperaturas do interior: valor que resume o resultado.
    double temperatureSum() const {
        double sum = 0.0;
        for (int i = 1; i <= N; ++i) {
            const double* row = T + idx(i, 0);
            for (int j = 1; j <= N; ++j) sum += row[j];
        }
        return sum;
    }

    void printField() const {
        std::printf("-----------Solução-----------\n");
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) std::printf("%.2f\t", T[idx(i, j)]);
            std::printf("\n");
        }
    }

    void writeVTK() const {
        exportToVTK(T, N, W, dx, dy, "temperatureField.vtk");

        std::vector<double> kField(cells());
        for (long p = 0; p < cells(); ++p) kField[p] = conductivity[material[p]];
        exportToVTK(kField.data(), N, W, dx, dy, "K_conductivity.vtk");
    }

private:
    const int N;     // células do interior por direção
    const int W;     // posições por linha no array (N + 2, com as células fantasma)
    const double dx, dy;
    const double tol;
    const long maxIter;

    // Todos os campos são arrays planos de W*W posições, linha por linha.
    std::vector<unsigned char> material;  // número do material (1..4) de cada célula
    std::vector<double> cE, cW, cN, cS;   // coeficientes já divididos por a_P
    std::vector<double> robinA, robinB;   // borda direita: T[i][N+1] = robinA[i]*T[i][N] + robinB[i]
    std::vector<double> bufferA, bufferB; // os dois campos de temperatura
    double* T;     // campo novo (escrito nesta iteração)
    double* Told;  // campo antigo (lido nesta iteração)
    double avgChange = 0.0;

    long cells() const { return static_cast<long>(N + 2) * (N + 2); }
    long idx(int i, int j) const { return static_cast<long>(i) * W + j; }

    // Quatro quadrantes, cada um com um material.
    void buildMaterialMap() {
        for (int i = 0; i < W; ++i) {
            for (int j = 0; j < W; ++j) {
                unsigned char m;
                if      (i <= N / 2 && j <= N / 2) m = 1;
                else if (i <= N / 2)               m = 2;
                else if (j <= N / 2)               m = 3;
                else                               m = 4;
                material[idx(i, j)] = m;
            }
        }
    }

    // Condutividade na face entre duas células: média harmônica se os materiais forem diferentes.
    static double faceConductivity(int mCenter, int mNeighbor) {
        double kc = conductivity[mCenter];
        double kn = conductivity[mNeighbor];
        return (mCenter != mNeighbor) ? 2.0 * kc * kn / (kc + kn) : kc;
    }

    // Pré-calcula os coeficientes do estêncil. Dividir por a_P aqui evita
    // uma divisão por célula a cada iteração e reduz de 5 para 4 os arrays lidos.
    void buildCoefficients() {
        for (int i = 1; i <= N; ++i) {
            for (int j = 1; j <= N; ++j) {
                int m = material[idx(i, j)];
                double aE = faceConductivity(m, material[idx(i, j + 1)]) * dy / dx;
                double aW = faceConductivity(m, material[idx(i, j - 1)]) * dy / dx;
                double aN = faceConductivity(m, material[idx(i - 1, j)]) * dx / dy;
                double aS = faceConductivity(m, material[idx(i + 1, j)]) * dx / dy;
                double aP = aE + aW + aN + aS;
                cE[idx(i, j)] = aE / aP;
                cW[idx(i, j)] = aW / aP;
                cN[idx(i, j)] = aN / aP;
                cS[idx(i, j)] = aS / aP;
            }
        }

        // Borda direita (Robin): (k*T[i][N] + h*dx*Tf) / (k + h*dx), separado em dois coeficientes.
        for (int i = 1; i <= N; ++i) {
            double k = conductivity[material[idx(i, N)]];
            robinA[i] = k / (k + h * dx);
            robinB[i] = h * dx * Tf / (k + h * dx);
        }
    }

    // Bordas superior e inferior (Dirichlet) são constantes: escritas uma vez só.
    void setDirichletRows(double* t) {
        for (int j = 0; j < W; ++j) {
            t[idx(0, j)]     = T_top;
            t[idx(N + 1, j)] = T_bottom;
        }
    }

    void applyBoundaryConditions(double* t) {
        for (int i = 1; i <= N; ++i) {
            double* row = t + idx(i, 0);
            row[0]     = row[1];                            // esquerda: Neumann (adiabática)
            row[N + 1] = robinA[i] * row[N] + robinB[i];    // direita: Robin (convecção)
        }
    }

    // Uma varredura de Jacobi no interior. A soma de |novo - antigo| é feita
    // no mesmo laço, sem um array de resíduo separado.
    double jacobiSweep() {
        double changeSum = 0.0;
        for (int i = 1; i <= N; ++i) {
            const double* up     = Told + idx(i - 1, 0);
            const double* center = Told + idx(i, 0);
            const double* down   = Told + idx(i + 1, 0);
            const double* e = &cE[idx(i, 0)];
            const double* w = &cW[idx(i, 0)];
            const double* n = &cN[idx(i, 0)];
            const double* s = &cS[idx(i, 0)];
            double* out = T + idx(i, 0);

            for (int j = 1; j <= N; ++j) {
                double v = e[j] * center[j + 1] + w[j] * center[j - 1]
                         + n[j] * up[j]         + s[j] * down[j];
                changeSum += std::fabs(v - center[j]);
                out[j] = v;
            }
        }
        return changeSum;
    }
};

static Settings parseArgs(int argc, char** argv) {
    Settings s;
    int positional = 0;
    for (int a = 1; a < argc; ++a) {
        if      (std::strcmp(argv[a], "--vtk") == 0)   s.vtk = true;
        else if (std::strcmp(argv[a], "--print") == 0) s.print = true;
        else if (positional == 0) { s.N = std::atoi(argv[a]);        ++positional; }
        else if (positional == 1) { s.tol = std::atof(argv[a]);      ++positional; }
        else if (positional == 2) { s.maxIter = std::atol(argv[a]);  ++positional; }
        else {
            std::fprintf(stderr, "Argumento desconhecido: %s\n", argv[a]);
            std::exit(1);
        }
    }
    if (s.N < 2 || s.tol <= 0.0 || s.maxIter < 1) {
        std::fprintf(stderr, "Uso: %s [N>=2] [tol>0] [maxIter>=1] [--vtk] [--print]\n", argv[0]);
        std::exit(1);
    }
    return s;
}

int main(int argc, char** argv) {
    using clock = std::chrono::steady_clock;
    auto seconds = [](clock::time_point a, clock::time_point b) {
        return std::chrono::duration<double>(b - a).count();
    };

    Settings settings = parseArgs(argc, argv);

    auto t0 = clock::now();
    HeatConductionSolver solver(settings);       // alocação, mapa de materiais e coeficientes
    auto t1 = clock::now();
    long iterations = solver.solve();            // laço iterativo
    auto t2 = clock::now();
    double sum = solver.temperatureSum();
    if (settings.print) solver.printField();
    if (settings.vtk)   solver.writeVTK();
    auto t3 = clock::now();

    double tInit = seconds(t0, t1), tSolve = seconds(t1, t2), tOut = seconds(t2, t3);
    double tTotal = seconds(t0, t3);

    std::printf("N               : %d\n", settings.N);
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
    return 0;
}
