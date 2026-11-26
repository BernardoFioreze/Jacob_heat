CXX      ?= g++
MPICXX   ?= mpicxx
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra
BIN      ?= heat
MPI_BIN  ?= heat_mpi

# Versão sequencial (alvo padrão: make)
$(BIN): main.cpp vtkWriter.cpp solverParameters.h vtkWriter.h
	$(CXX) $(CXXFLAGS) -o $@ main.cpp vtkWriter.cpp

# Versão paralela com MPI: make heat_mpi
$(MPI_BIN): heat_mpi.cpp vtkWriter.cpp solverParameters.h vtkWriter.h
	$(MPICXX) $(CXXFLAGS) -o $@ heat_mpi.cpp vtkWriter.cpp

clean:
	rm -f heat heat_mpi *.vtk

.PHONY: clean
