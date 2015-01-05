#include "Bspline.h"

int main() {
    vector<double> X = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0 };

    Bspline bspline( X, 3, 1 );

    bspline.bspline();    
}
