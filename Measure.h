#ifndef __MEASURE__
#define __MEASURE__

#include "Profile.h"
#include "Outcome.h"
#include <cmath>
#include <cstring>
#include <algorithm>
using namespace std;

inline double entropy( int up, int dn ) { 
	if ( up == 0 ) return 0.0;
	else return - double(up) / dn * log2( double(up) / dn );
} 

double chiSquare(	vector<int> A, const int numTypeA, 
					vector<int> B, const int numTypeB,
					Outcome outcome ) {

	const int N = (int)outcome.size();
	const int maxState = numTypeA * numTypeB;

	int freq[ outcome.getNumTypes() ][ maxState ];
	int rowsum[ outcome.getNumTypes() ];
	int colsum[ maxState ];

	memset( freq, 0, sizeof freq );
	memset(rowsum, 0, sizeof rowsum);
	memset(colsum, 0, sizeof colsum);

	for( size_t i = 0 ; i < N ; ++i ) {
		freq[ outcome[i] ][ A[i] * numTypeB + B[i] ]++;	
		rowsum[ outcome[i] ]++;
		colsum[ A[i] * numTypeB + B[i] ]++;
	}

	double ret = 0.0;


	for( int i = 0 ; i < outcome.getNumTypes() ; ++i ) {
		for( int j = 0 ; j < maxState ; ++j ) {
			double expected = rowsum[i] * colsum[j] / (double)N;
			if( fabs(expected) < 1e-9 ) continue;
			ret += (expected-freq[i][j]) * (expected-freq[i][j]) / expected; 
		}
	}
	return ret;
	
}


double mutualInformation(   vector<int> A, const int numTypeA, 
                            vector<int> B, const int numTypeB, 
                            Outcome outcome ) {

	const int maxState = numTypeA * numTypeB;

	int freq[ outcome.getNumTypes() ][ maxState ];

	memset( freq, 0, sizeof freq );

	for( size_t i = 0 ; i < outcome.size() ; ++i ) {
		freq[ outcome[i] ][ A[i] * numTypeB + B[i] ]++;	
	}

	double H_Y = 0;
    for( int i = 0 ; i < outcome.getNumTypes() ; ++i ) {
        H_Y += entropy( outcome.getNumSubjects(i), outcome.size() );
    }

	double H_X = 0;
	for( int i = 0 ; i < maxState ; ++i ) {
		int colsum = 0;
		for( int j = 0 ; j < outcome.getNumTypes() ; ++j ) {
			colsum += freq[j][i];
		}
		H_X += entropy( colsum, outcome.size() );
	}

	double H_XY = 0;
	for( int i = 0 ; i < outcome.getNumTypes() ; ++i ) {
		for( int j = 0 ; j < maxState ; ++j ) {
			H_XY += entropy( freq[i][j], outcome.size() );
		}
	}
	
	return ( H_X + H_Y - H_XY );
	return 0;
}

#endif
