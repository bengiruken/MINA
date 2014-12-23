#ifndef __NETWORK__
#define __NETWORK__
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
#include "Outcome.h"
#include "Profile.h"
#include "Measure.h"
#include "Param.h"
using namespace std;

double getOutcomeAssociationNetworkThreshold( Profile &profile, Outcome &outcome, const int numPermute ) {
    srand( time(NULL) );
    vector<Outcome> outcomes( numPermute, outcome );
    vector<double> sumMI( numPermute, 0 );
    for( size_t  i = 0 ; i < numPermute ; ++i ) {
        outcomes[i].permute();
    }

    double maxi = 0;

    const long long totalIteration = (long long) profile.getNumFeatures() * 
        (profile.getNumFeatures()-1) / 2;

    showProgress( 0, 1, true );
    long long iteration = 0;
    for( size_t i = 0 ; i < profile.getNumFeatures() ; ++i ) {
        for( size_t j = i+1 ; j < profile.getNumFeatures() ; ++j ) {

            # pragma omp for
            for( int iter = 0 ; iter < numPermute ; ++iter ) {
                sumMI[iter] = getOutcomeAssociation(    
                    profile[i], profile.getNumTypes(), 
                    profile[j], profile.getNumTypes(),
                    outcomes[iter] );
            }
            showProgress( ++iteration, totalIteration );
            double avg = accumulate( sumMI.begin(), sumMI.end(), 0.0 ) / numPermute;
            maxi = max( maxi, avg );
        }
    }

    return maxi;
}

double getInteractionNetworkThreshold( Profile &profile, const int numPermute ) {
    srand( time(NULL) );
    vector<double> sumMI( numPermute, 0 );

    double maxi = 0;

    const long long totalIteration = (long long) profile.getNumFeatures() * 
        (profile.getNumFeatures()-1) / 2;

    showProgress( 0, 1, true );
    long long iteration = 0;
    for( size_t i = 0 ; i < profile.getNumFeatures() ; ++i ) {
        vector< vector<int> > pi = vector< vector< int > >( numPermute, profile[i] );
        for( int iter = 0 ; iter < numPermute ; ++iter ) {
            random_shuffle( pi[iter].begin(), pi[iter].end() );
        }

        for( size_t j = i+1 ; j < profile.getNumFeatures() ; ++j ) {

            # pragma omp for
            for( int iter = 0 ; iter < numPermute ; ++iter ) {
                sumMI[iter] = getInteraction( 
                    pi[iter], profile.getNumTypes(),
                    profile[j], profile.getNumTypes() );
            }
            showProgress( ++iteration, totalIteration );
            double avg = accumulate( sumMI.begin(), sumMI.end(), 0.0 ) / numPermute;

            if( maxi < avg ) {
                cerr << "[DEBUG]" << "maxi : " << maxi << " at (";
                cerr << i << "," << j << ")" << endl;
            }
            maxi = max( maxi, avg );
        }
    }

    return maxi;
}

typedef vector<int> EdgeList;
typedef vector< EdgeList > Network;

Network getAssocationNetwork( Profile &profile, Outcome &outcome, double thresholds, double alpha ) {
    const double strictThreshold = (1+alpha) * thresholds;
    const int numFeatures = profile.getNumFeatures();
    const long long totalIteration = (long long) numFeatures * (numFeatures-1) / 2;

    long long iteration = 0;

    showProgress( 0, totalIteration, true );

    Network network( numFeatures, EdgeList(0) );

    for( int i = 0 ; i < numFeatures ; ++i ) {
        for( int j = i+1 ; j < numFeatures ; ++j ) {
            double value = getOutcomeAssociation( 
                profile[i], profile.getNumTypes(), 
                profile[j], profile.getNumTypes(), outcome );
            if( value > strictThreshold ) {
                network[i].push_back(j);
            }
        }
    }
    return network;
}


Network getInteractionNetwork( Profile &profile, double thresholds, double alpha ) {
    const double strictThreshold = (1+alpha) * thresholds;
    const int numFeatures = profile.getNumFeatures();
    const long long totalIteration = (long long) numFeatures * (numFeatures-1) / 2;

    long long iteration = 0;

    showProgress( 0, totalIteration, true );

    Network network( numFeatures, EdgeList(0) );

    for( int i = 0 ; i < numFeatures ; ++i ) {
        for( int j = i+1 ; j < numFeatures ; ++j ) {
            double value = getInteraction(
                profile[i], profile.getNumTypes(), 
                profile[j], profile.getNumTypes() );
            if( value > strictThreshold ) {
                network[i].push_back(j);
            }
        }
    }
    return network;
}


#endif
