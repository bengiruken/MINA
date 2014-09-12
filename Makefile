all: main.cc
	g++ -fopenmp -O3 main.cc -o main_th
