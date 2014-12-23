all: main.cc
	g++ -fopenmp -O3 main.cc -o bin/MINA

branch: main.cc
	g++ -fopenmp -O3 test.cc -o bin/PIONA
