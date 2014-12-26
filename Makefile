all: mina.cc piona.cc
	g++ -fopenmp -O3 mina.cc -o bin/MINA
	g++ -fopenmp -O3 piona.cc -o bin/PIONA
