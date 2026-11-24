#include "vtkWriter.h"
#include <cstdio>

void exportToVTK(const double* field, int N, int stride, double dx, double dy,
                 const std::string& filename) {
    FILE* f = std::fopen(filename.c_str(), "w");
    if (!f) {
        std::fprintf(stderr, "Erro: não foi possível abrir %s para escrita.\n", filename.c_str());
        return;
    }

    std::fprintf(f, "# vtk DataFile Version 3.0\n");
    std::fprintf(f, "Generated from Composite Slab Solver\n");
    std::fprintf(f, "ASCII\n");
    std::fprintf(f, "DATASET STRUCTURED_POINTS\n");
    std::fprintf(f, "DIMENSIONS %d %d 1\n", N, N);
    std::fprintf(f, "ORIGIN 0 0 0\n");
    std::fprintf(f, "SPACING %g %g 1\n", dx, dy);
    std::fprintf(f, "POINT_DATA %ld\n", static_cast<long>(N) * N);
    std::fprintf(f, "SCALARS field double\n");
    std::fprintf(f, "LOOKUP_TABLE default\n");

    for (int i = 1; i <= N; ++i) {
        const double* row = field + static_cast<long>(i) * stride;
        for (int j = 1; j <= N; ++j) {
            std::fprintf(f, "%.3f ", row[j]);
        }
        std::fputc('\n', f);
    }

    std::fclose(f);
    std::printf("Arquivo escrito: %s\n", filename.c_str());
}
