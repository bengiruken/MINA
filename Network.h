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
            maxi = max( maxi, avg );
        }
    }

    return maxi;
}

typedef pair<int,double> Edge;
typedef vector<Edge> EdgeList;
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
                network[i].push_back(Edge(j,value));
            }
            showProgress(++iteration, totalIteration);
        }
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
                network[i].push_back(Edge(j,value));
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
            out << " " << here[j].first << " " << here[j].second;
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
        double weight;
        EdgeList edgeList;
        while( sin >> vertex >> weight ) {
            edgeList.push_back(Edge(vertex,weight));
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
        EdgeList::const_iterator it = mine[i].begin();
        EdgeList::const_iterator that = other[i].begin();

        while( it != mine[i].end() && that != other[i].end() ) {
            if( *it == *that ) {
                ++it;
                ++that;
            }
            else if( *it < *that ) {
                network[i].push_back(*it);
                ++it;
            }
            else {
                ++that;
            }
                
        }

    }
    return network;
}

Network filterNetwork( const Network &network, const double &cutoff ) {
    const int N = network.size();

    Network ret(N,EdgeList(0));
  
    for( int i = 0 ; i < N ; ++i ) {
        for( auto it = network[i].begin() ; it != network[i].end() ; ++it ) {
            if( it->second > cutoff ) ret[i].push_back(*it);
        }
    }

    return ret;
}

#endif
