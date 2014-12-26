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

Network getAssocationNetwork( Profile &profile, Outcome &outcome, double thresholds, double alpha = 0.0 ) {
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
            showProgress(++iteration, totalIteration);
        }
        // fprintf( stderr, "[DEBUG] network[%d] = %d\n", i, network[i].size() );
    }
    return network;
}


Network getInteractionNetwork( Profile &profile, double thresholds, double alpha = 0.0 ) {
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
            showProgress(++iteration, totalIteration);
        }
    }
    return network;
}

void saveNetwork( const Network &network, const char *fileName ) {
    ofstream out(fileName);
    const int N = network.size();

    int numEdges = 0;
    for( int i = 0 ; i < N ; ++i ) {
        const EdgeList &here = network[i];
        const int M = here.size();
        out << i << ":";
        for( int j = 0 ; j < M ; ++j ) {
            out << " " << here[j];
        }         
        numEdges += M;
        out << endl;
    }
    out.close();

    cout << "[SAVE NETWORK]" << endl;
    cout << "File name = " << fileName << endl;
    cout << "Num vertices = " << N << endl;
    cout << "Num edges = " << numEdges << endl;

}

Network loadNetwork( const char *fileName ) {
    Network network;
    string line;
    ifstream fin(fileName);
   
    int numVertices = 0;
    int numEdges = 0;
    while( getline( fin, line ) ) {
        ++numVertices;
        istringstream sin( line );

        string tok;

        sin >> tok;

        int vertex;

        EdgeList edgeList;
        while( sin >> vertex ) {
            edgeList.push_back(vertex);
        }
    
        numEdges += edgeList.size();
        network.push_back(edgeList);
    }
    cout << "[LOAD NETWORK]" << endl;
    cout << "File name = " << fileName << endl;
    cout << "Num vertices = " << numVertices << endl;
    cout << "Num edges = " << numEdges << endl;
    fin.close();

    return network;
}


Network getDifferenceNetwork( const Network &mine, const Network &other ) {
    assert( mine.size() == other.size() );
    Network network( mine.size(), EdgeList(0) );
    const int N = mine.size();

    for( int i = 0 ; i < N ; ++i ) {
        EdgeList edgeList( mine[i].size() + other[i].size(), 0 );
        EdgeList::iterator it = set_difference( mine[i].begin(), mine[i].end(), 
            other[i].begin(), other[i].end(), edgeList.begin() );
        edgeList.resize(it-edgeList.begin());

        network[i] = edgeList;
    }
    return network;
}

#endif
