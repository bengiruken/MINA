all: mina.cc piona.cc
	g++ -fopenmp -O3 --std=c++11 mina.cc -o bin/MINA
	g++ -fopenmp -O3 --std=c++11 piona.cc -o bin/PIONA
