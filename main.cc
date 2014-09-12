#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <numeric>
#include "omp.h"
#include "Utility.h"
#include "Outcome.h"
#include "Profile.h"
#include "Measure.h"

using namespace std;

template <typename T>
void permutationTest( Profile<T> &profile, Outcome &outcome, vector<string> &genenames, FILE *fp, const int numPermute = int(1e4) ) {
	double mival;
	for( size_t i = 0 ; i < profile.getNumFeatures() ; ++i ) {
		cerr << i+1 << "/" << profile.getNumFeatures() << endl;
		for( size_t j = 0 ; j < profile.getNumFeatures() ; ++j ) {
			fscanf(fp,"%lf", &mival);
			if( mival < 0.1 ||  i > j ) continue;

			vector<int> numover( numPermute, 0 );
			vector<Outcome> outcomes( numPermute, outcome );

			# pragma omp for
			for( int iter = 0 ; iter < numPermute ; ++iter ) {
				outcomes[iter].permute();
				numover[iter] = mival <= mutualInformation( profile[i], profile.getNumTypes(), profile[j], profile.getNumTypes(), outcomes[iter] ); 
			}

			int sumnumover = accumulate( numover.begin(), numover.end(), 0 );
			printf("%s\t%s\t%f\t%f\n", genenames[i].c_str(), genenames[j].c_str(), mival, sumnumover / double(numPermute) );
		}
	}

}

template <typename T>
void getThreshold( Profile<T> &profile, Outcome &outcome, const int numPermute, vector<bool> &selected ) {
	srand( time(NULL) );
	vector<Outcome> outcomes( numPermute, outcome );
	vector<double> sumMI( numPermute, 0 );
	for( size_t  i = 0 ; i < numPermute ; ++i ) {
		outcomes[i].permute(selected);
	}
	// int total = 0;
	
	double maxi = -1e100;
	for( size_t i = 0 ; i < profile.getNumFeatures() ; ++i ) {
		for( size_t j = i+1 ; j < profile.getNumFeatures() ; ++j ) {

			# pragma omp for
			for( int iter = 0 ; iter < numPermute ; ++iter ) {
				sumMI[iter] = mutualInformation( profile[i], profile.getNumTypes(), profile[j], profile.getNumTypes(), outcomes[iter], selected );
			}
			double avg = accumulate( sumMI.begin(), sumMI.end(), 0.0 ) / numPermute;
			// printf("%.4f ", avg );
			maxi = max( maxi, avg );
		}
		//printf("\n");
		// total += profile.getNumFeatures() - (i+1);
		// cerr << "max : " << *max_element( sumMI.begin(), sumMI.end() ) / total << endl;
		
		if( i % 100 == 0 ) {
			cerr << i+1 << "/" << profile.getNumFeatures() << endl;
			cerr << "max : " << maxi << endl;
		}
	}
	cout << maxi << " ";
	
	/*
	for( int i = 0 ; i < numPermute ; ++i ) {
		cout << i << "\t" << sumMI[i] / total << endl;
	}*/
}


template <typename T>
void getCompleteGraph( Profile<T> &profile, Outcome &outcome ) {
	outcome.permute();

	for( size_t i = 0 ; i < profile.getNumFeatures() ; ++i ) {
		for( size_t j = 0 ; j < profile.getNumFeatures() ; ++j ) {
			if( i == j ) {
				printf("0.0000\t");
				continue;
			}
			else {
				double mival = mutualInformation( profile[i], profile.getNumTypes(), profile[j], profile.getNumTypes(), outcome );
				printf("%.4f\t", mival);
			}
		}
		printf("\n");
	}
}


template <typename T>
void getCompleteGraphForCV( Profile<T> &profile, Outcome &outcome, vector<bool> &selected ) {

	for( size_t i = 0 ; i < profile.getNumFeatures() ; ++i ) {
		for( size_t j = 0 ; j < profile.getNumFeatures() ; ++j ) {
			if( i == j ) {
				printf("0\t");
				continue;
			}
			else {
				double mival = mutualInformation( profile[i], profile.getNumTypes(), profile[j], profile.getNumTypes(), outcome, selected );
				printf("%.4f\t", mival);
			}
		}
		printf("\n");
	}
}


int main(int argv, char *argc[]) {
	Outcome outcome("../data20140108/outcome.dat");
	vector<string> genenames = readSymbols("../data20140108/genename.txt");
	/*
	if( strcmp(argc[1],"int") == 0 ) {
		Profile<int> profile(argc[2]);
		// permutationTest( profile, outcome, genenames, fp );
		profile.toDiscrete( genenames );
	}
	else {
		Profile<double> profile(argc[2]);
		// permutationTest( profile, outcome, genenames, fp );
		profile.toDiscrete( genenames );
	}*/
	
	
	/*
	FILE *fp = fopen(argc[3],"r");
	fclose(fp);
	if( strcmp(argc[1],"int") == 0 ) {
		Profile<int> profile(argc[2]);
		permutationTest( profile, outcome, genenames, fp );
	}
	else {
		Profile<double> profile(argc[2]);
		permutationTest( profile, outcome, genenames, fp );
	} */

	FILE *fp = 	fopen("../foldIdx.txt","r");

	int idx = 0;

	//sscanf(argc[3],"%d", &idx);

	vector<bool> selected;
	int ino = 0;
	while( fscanf(fp,"%d", &ino) == 1 ) {
		selected.push_back( ino != idx );	
	}
	/*
	if( strcmp(argc[1],"int") == 0 ) {
		Profile<int> profile(argc[2]);
		getCompleteGraphForCV( profile, outcome, selected );
	}
	else {
		Profile<double> profile(argc[2]);
		getCompleteGraphForCV( profile, outcome, selected );
	}*/
	if( strcmp(argc[1],"int") == 0 ) {
		Profile<int> profile(argc[2]);
		getThreshold( profile, outcome, 30, selected );
	}
	else {
		Profile<double> profile(argc[2]);
		getThreshold( profile, outcome, 30, selected );
	}
	cerr << argc[2] << "is done!" << endl;
	return 0;
}

