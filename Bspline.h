#ifndef __BSPLINE__
#define __BSPLINE__
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;

struct Bspline {

    vector<double> X;
    int M, K;
    int N;
    Bspline( vector<double> _X, int _M, int _K ) : X(_X), M(_M), K(_K) {
        N = (int)X.size();
    }

    double T( const int i ) {
        if( i < K ) return 0.0;
        else if( i < M ) return i - K + 1.0;
        else return M - K + 1.0;
    }

    vector< double > normalize() {
        const int N = X.size();
        vector< double > Z( N, 0 );

        auto minX = *min_element( X.begin(), X.end() );
        auto maxX = *max_element( X.begin(), X.end() );

        double disp = (M-K+1) / (maxX - minX);
        for( int i = 0 ; i < N ; ++i ) {
            Z[i] = (X[i] - minX) * disp;
        }
        return Z;
    }

    double cache[100][100];
    double B( const double Z, const int i, const int k ) {
        double &ret = cache[i][k];
        if( k == 1 ) {
            if( T(i) <= Z && Z < T(i+1) ) ret = 1.0;
            else if( fabs( Z-T(i) ) < 1e-10 && i+1 >= M ) ret = 1.0;
            else ret = 0.0;

            return ret;
        }
        else {
            if(ret >= 0.0) return ret;
            ret = 0.0;

            double d1 = T(i+k-1) - T(i);
            double d2 = T(i+k) - T(i+1);
            if( fabs(d1) > 1e-10 ) 
                ret += B(Z,i,k-1) * (Z-T(i))/ d1;
            if( fabs(d2) > 1e-10 ) 
                ret += B(Z,i+1,k-1) * (T(i+k)-Z) / d2;
            return ret;
        }
    }

    vector< double > getBX( const double Z ) {
        for( int i = 0 ; i < 100 ; ++i ) {
            for( int j = 0 ; j < 100 ; ++j ) {
                cache[i][j] = -1e100;
            }
        }

        vector<double> recent(M+1,0);

        for( int i = 0 ; i <= M ; ++i ) {
            recent[i] = B(Z,i,K);
        }

        for( int k = 1 ; k <= 2 ; ++k ) {
            for( int i = 0 ; i <= M ; ++i ) {
                cerr << cache[i][k] << " ";
            }
            cerr << endl;
        }
        return recent;
    }

    vector< vector<double> > bspline() {
        vector< vector<double> > BX( N, vector<double>(M,0) );

        auto BZ = normalize();

        for( int i = 0 ; i < N ; ++i ) {
            BX[i] = getBX( BZ[i] );
        }
        return BX;
    }
};

#endif
