all: main.cc
	g++ -fopenmp -O3 main.cc -o bin/MINA
