CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++11 -Wall -Wextra

poisson: src/poisson.cpp src/vtk_writer.hpp
	$(CXX) $(CXXFLAGS) -o poisson src/poisson.cpp

rodar: poisson
	./poisson --dim 3 --n 49 --caso seno --saida saida

limpar:
	rm -f poisson

.PHONY: rodar limpar
