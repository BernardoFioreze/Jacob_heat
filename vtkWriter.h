#pragma once
#include <string>

/*
 * Exporta um campo escalar 2D para o formato VTK legado (abre no ParaView).
 *
 * O campo é um array plano, linha por linha, com (N+2) x (N+2) posições:
 * o interior N x N mais uma camada de células fantasma em volta.
 * Só o interior é escrito no arquivo.
 *
 * field    : início do array
 * N        : número de células do interior em cada direção
 * stride   : número de posições por linha no array (N + 2)
 * dx, dy   : espaçamento da grade
 * filename : nome do arquivo .vtk de saída
 */
void exportToVTK(const double* field, int N, int stride, double dx, double dy,
                 const std::string& filename);
