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
double mutualInformation( vector<int> A, const int numTypeA, vector<int> B, const int numTypeB, Outcome outcome, vector<bool> selected ) {

	const int maxState = numTypeA * numTypeB;

	const int TOT = count(selected.begin(),selected.end(),true);
	int freq[ outcome.getNumTypes() ][ maxState ];

	memset( freq, 0, sizeof freq );

	for( size_t i = 0 ; i < outcome.size() ; ++i ) if( selected[i] ) {
		freq[ outcome[i] ][ A[i] * numTypeA + B[i] ]++;	
	}
	double H_Y = entropy( outcome.getNumSubjects(0) , TOT ) + entropy( outcome.getNumSubjects(1), TOT );

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

	return ( H_X + H_Y - H_XY ) / H_Y;
	return 0;
}



double mutualInformation( vector<int> A, const int numTypeA, vector<int> B, const int numTypeB, Outcome outcome) {

	const int maxState = numTypeA * numTypeB;

	int freq[ outcome.getNumTypes() ][ maxState ];

	memset( freq, 0, sizeof freq );

	for( size_t i = 0 ; i < outcome.size() ; ++i ) {
		freq[ outcome[i] ][ A[i] * numTypeA + B[i] ]++;	
	}
	double H_Y = entropy( outcome.getNumSubjects(0) , outcome.size() ) + entropy( outcome.getNumSubjects(1), outcome.size() );

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

	return ( H_X + H_Y - H_XY ) / H_Y;
	return 0;
}

double permutationTest( vector<int> &A, const int numTypeA, vector<int> &B, const int numTypeB, Outcome outcome, int numIteration ) {
	double realValue = mutualInformation( A, numTypeA, B, numTypeB, outcome );	

	int numFailure;
	for( int i = 0 ; i < numIteration ; ++i ) {
		outcome.permute();
		double randomValue = mutualInformation( A, numTypeA, B, numTypeB, outcome );
		if( randomValue >= realValue ) ++numFailure;
	}

	return numFailure / double(numIteration);
}

#endif
