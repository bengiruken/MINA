#ifndef __MEASURE__
#define __MEASURE__

#include "Profile.h"
#include "Outcome.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <iomanip>

using namespace std;

inline double plogp( int up, int dn ) { 
	if ( up == 0 ) return 0.0;
	else return double(up) / dn * log2( double(up) / dn );
} 


double getInteraction(  vector<int> A, const int numTypeA, 
                        vector<int> B, const int numTypeB ) {

    const int N = (int)A.size();

	vector< int > freq(numTypeA*numTypeB);
	vector< int > freqX(numTypeA, 0);
	vector< int > freqY(numTypeB, 0);

	for( int i = 0 ; i < N ; ++i ) {
        freq[ A[i] * numTypeB + B[i] ]++;
		freqX[A[i]]++;
		freqY[B[i]]++;
	}

    double H_X = 0;
	for (auto x : freqX) {
		H_X += -plogp(x, N);
	}

	double H_Y = 0;
	for (auto y : freqY) {
		H_Y += -plogp(y, N);
	}

	
    double H_XY = 0;
	for (auto xy : freq) {
		H_XY += -plogp(xy, N);
	}
	
	return (H_X + H_Y - H_XY) / min(H_X,H_Y);
}

double getOutcomeAssociation(   vector<int> A, const int numTypeA, 
                                vector<int> B, const int numTypeB, 
                                Outcome outcome ) {

	const int maxState = numTypeA * numTypeB;

	vector< vector<int> > freq( outcome.getNumTypes(), vector<int>( maxState, 0));
	vector< int > freqCombination(maxState, 0);
	vector< int > freqOutcome(outcome.getNumTypes(), 0);
	for( size_t i = 0 ; i < outcome.size() ; ++i ) {
		freq[ outcome[i] ][ A[i] * numTypeB + B[i] ]++;
		freqCombination[A[i] * numTypeB + B[i]]++;
		freqOutcome[outcome[i]]++;
	}

	double H_X = 0;
	for ( auto x : freqOutcome ) {
        H_X += -plogp( x, outcome.size() );
    }

	double H_Y = 0;
	for (auto x : freqCombination) {
		H_Y += -plogp(x, outcome.size());
	}
	double H_XY = 0;
	for (auto r : freq) {
		for (auto x : r) {
			H_XY += -plogp(x, outcome.size());
		}
	}

	return (H_X + H_Y - H_XY) / min(H_X,H_Y);


	return 0;
}

double entropy( const vector<int> )
#endif
