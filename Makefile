all: main.cc
	g++ -std=c++11 -fopenmp -O3 main.cc -o bin/MINA
