#pragma once

// Parâmetros físicos fixos do problema.
// O tamanho da grade (N), a tolerância e o limite de iterações são lidos
// da linha de comando (ver main.cpp), com os valores padrão abaixo.

// Geometria da placa (m)
const double xLength = 0.1;
const double yLength = 0.1;

// Condutividade térmica de cada material (W/m.K), indexada pelo número do material (1..4).
// O índice 0 não é usado.
const double conductivity[5] = {0.0, 50.0, 100.0, 150.0, 200.0};

// Condições de contorno
const double T_top    = 100.0;  // Dirichlet na borda superior (°C)
const double T_bottom = 200.0;  // Dirichlet na borda inferior (°C)
const double h  = 50.0;         // coeficiente de convecção da borda direita (W/m².K)
const double Tf = 50.0;         // temperatura do fluido na borda direita (°C)
const double T_initial = 1.0;   // chute inicial do campo de temperatura

// Valores padrão da linha de comando
const int    default_N       = 16;
const double default_tol     = 1e-6;
const long   default_maxIter = 100000000;  // só uma trava de segurança: a parada normal é pela tolerância
