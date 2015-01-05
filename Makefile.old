all: mina.cc piona.cc
	g++ -fopenmp -O3 --std=c++11 mina.cc -o bin/MINA
	g++ -fopenmp -O3 --std=c++11 piona.cc -o bin/PIONA
	g++ -fopenmp -O3 --std=c++11 toy_example.cc -o bin/toy_example

bspline: bspline.cc
	g++ -fopenmp -O3 --std=c++11 bspline.cc -o bin/bspline
