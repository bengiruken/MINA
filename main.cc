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
#include "Param.h"
using namespace std;

double getThreshold(    Profile &profile, Outcome &outcome, 
                        const int numPermute ) {
    srand( time(NULL) );
    vector<Outcome> outcomes( numPermute, outcome );
    vector<double> sumMI( numPermute, 0 );
    for( size_t  i = 0 ; i < numPermute ; ++i ) {
        outcomes[i].permute();
    }

    double maxi = 0;

    const long long totalIteration = (long long) profile.getNumFeatures() * 
                                                 (profile.getNumFeatures()-1) / 2;

    long long iteration = 0;
    for( size_t i = 0 ; i < profile.getNumFeatures() ; ++i ) {
        for( size_t j = i+1 ; j < profile.getNumFeatures() ; ++j ) {

            # pragma omp for
            for( int iter = 0 ; iter < numPermute ; ++iter ) {
                sumMI[iter] = mutualInformation(    profile[i], profile.getNumTypes(), 
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

void getEdges(  const vector<string> genenames, const Param &param, 
                vector<Profile> &profiles, Outcome &outcome, 
                const vector<double>  &thresholds ) {

    const size_t numProfiles = param.profiles.size();

    vector<string> fileName;
    vector<FILE *> fp;
    for( size_t i = 0 ; i < numProfiles ; ++i ) {
        fileName.push_back( "output/"+base_name(param.profiles[i]) );
        fp.push_back( fopen(fileName.back().c_str(),"w") );
    }

    FILE *outUni = fopen("output/union.txt","w");
    FILE *outInt = fopen("output/inter.txt","w");

    const vector<double> &alpha = param.alpha;

    const size_t numFeatures = profiles.front().getNumFeatures();

    const long long totalIteration = (long long) numFeatures * 
                                                 (numFeatures-1) / 2;
    long long iteration = 0;

    showProgress( 0, totalIteration, true );
    for( size_t i = 0 ; i < numFeatures ; ++i ) for( size_t j = i + 1 ; j < numFeatures ; ++j ) {
        vector<double> values( numProfiles, 0.0 );
        for( int ord = 0 ; ord < numProfiles ; ++ord ) {
            Profile &P = profiles[ord];
            values[ord] = mutualInformation(  P[i], P.getNumTypes(), 
                                              P[j], P.getNumTypes(),
                                              outcome );
        }

        bool isInter = true;
        bool isUnion = false;
        for( int ord = 0 ; ord < numProfiles ; ++ord ) {
            if( values[ord] <= alpha[1] * thresholds[ord] ) {
                isInter = false; 
            }
            
            if( values[ord] > alpha[2] * thresholds[ord] ) {
                isUnion = true;
            }
            if( values[ord] > alpha[1] * thresholds[ord] ) {
                fprintf( fp[ord], "%s\t%s\t%.4f\n",   genenames[i].c_str(),
                                                      genenames[j].c_str(),
                                                      values[ord] );
            }

        }

        if( isInter ) {
            fprintf( outInt, "%s\t%s", genenames[i].c_str(), genenames[j].c_str() );
            for( int ord = 0 ; ord < numProfiles ; ++ord ) {
                fprintf( outInt, "\t%.4f", values[ord] );
            }
            fprintf( outInt, "\n" );

        }
        
        if( isUnion ) {
            fprintf( outUni, "%s\t%s", genenames[i].c_str(), genenames[j].c_str() );
            for( int ord = 0 ; ord < numProfiles ; ++ord ) {
                if( values[ord] > alpha[2] * thresholds[ord] ) {
                    fprintf( outUni, "\t%.4f", values[ord] );
                }
                else {
                    fprintf( outUni, "\t0.0000", values[ord] );
                }
            }
            fprintf( outUni, "\n" );
        }


        showProgress( ++iteration, totalIteration );
    }

    fclose( outUni );
    fclose( outInt );
}

int main(int argv, char *argc[]) {

    cerr << "Starting MINA..." << endl;
    Param param;
    cerr << "Getting parameters...";
    ofstream outPa("output/param.txt");
    param.getParamInfo(outPa);
    outPa.close();
    cerr << "DONE!" << endl;

    Outcome outcome(param.clinical.c_str());
    vector<string> genenames = readSymbols(param.geneInfo.c_str());
    vector<Profile> profiles;
    vector<double> thresholds;


    cerr << "Reading profiles and get threshold value..." << endl;
    ofstream outTh("output/threshold.txt");

    for( size_t i = 0 ; i < param.profiles.size() ; ++i ) {
        profiles.push_back( Profile( param.profiles[i].c_str(), 5 ) );
        thresholds.push_back( getThreshold( profiles.back(), outcome, param.maxPerm ) );
        outTh << param.profiles[i] << "\t" << thresholds.back() << endl;
    }
    cerr << " DONE!" << endl;
    outTh.close();

    cerr << "Generating networks..." << endl;
    getEdges( genenames, param, profiles, outcome, thresholds );
    cerr <<  " DONE!" << endl;

    cerr << "Finish! You can check the results into output folder!" << endl;
    return 0;
}

