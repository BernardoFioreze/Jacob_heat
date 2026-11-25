CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra
BIN      ?= heat

$(BIN): main.cpp vtkWriter.cpp solverParameters.h vtkWriter.h
	$(CXX) $(CXXFLAGS) -o $@ main.cpp vtkWriter.cpp

clean:
	rm -f heat *.vtk

.PHONY: clean
