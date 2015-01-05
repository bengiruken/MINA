all: bspline.cc
	g++ -fopenmp -O3 --std=c++11 bspline.cc -o bin/bspline
	./bin/bspline
