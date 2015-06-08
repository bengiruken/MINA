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
#include "FrequencyTable.h"
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

    showProgress( 0, 1, true );
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
        fileName.push_back( param.outpath+base_name(param.profiles[i]) );
        fp.push_back( fopen(fileName.back().c_str(),"w") );
    }


    FILE *outUni = fopen((param.outpath+"union.txt").c_str(),"w");
    FILE *outInt = fopen((param.outpath+"output/inter.txt").c_str(),"w");

    double alpha = param.alpha;


    const size_t numFeatures = profiles.front().getNumFeatures();


    const long long totalIteration = (long long) numFeatures * 
                                                 (numFeatures-1) / 2;
    long long iteration = 0;

    showProgress( 0, totalIteration, true );
    for( size_t i = 0 ; i < numFeatures ; ++i ) for( size_t j = i + 1 ; j < numFeatures ; ++j ) {
        vector<double> values( numProfiles, 0.0 );

        # pragma omp for
        for( int ord = 0 ; ord < numProfiles ; ++ord ) {
            Profile &P = profiles[ord];
            values[ord] = mutualInformation(  P[i], P.getNumTypes(), 
                                              P[j], P.getNumTypes(),
                                              outcome );
        }

        bool isInter = true;
        bool isUnion = false;


        for( int ord = 0 ; ord < numProfiles ; ++ord ) {
            if( values[ord] <= alpha * thresholds[ord] ) {
                isInter = false; 
            }
            
            if( values[ord] > alpha * thresholds[ord] ) {
                isUnion = true;
            }
            if( values[ord] > alpha * thresholds[ord] ) {
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
                if( values[ord] > alpha * thresholds[ord] ) {
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

void getDist(  const Param &param, vector<Profile> &profiles, Outcome &outcome ) {
    const size_t numProfiles = param.profiles.size();

    vector<string> fileName;
    for( size_t i = 0 ; i < numProfiles ; ++i ) {
        fileName.push_back( param.outpath+base_name(param.profiles[i]) );
    }

    vector<FrequencyTable> ft( numProfiles );


    const size_t numFeatures = profiles.front().getNumFeatures();


    const long long totalIteration = (long long) numFeatures * 
                                                 (numFeatures-1) / 2;
    long long iteration = 0;

    showProgress( 0, totalIteration, true );
    for( size_t i = 0 ; i < numFeatures ; ++i ) for( size_t j = i + 1 ; j < numFeatures ; ++j ) {
        vector<double> values( numProfiles, 0.0 );

        # pragma omp for
        for( int ord = 0 ; ord < numProfiles ; ++ord ) {
            Profile &P = profiles[ord];
            values[ord] = mutualInformation(  P[i], P.getNumTypes(), 
                                              P[j], P.getNumTypes(),
                                              outcome );
            ft[ord].put(values[ord]);
        }



        showProgress( ++iteration, totalIteration );
    }

    for( int ord = 0 ; ord < numProfiles ; ++ord ) {
        ft[ord].output(fileName[ord].c_str(), 30, param.distLo, param.distHi );
    }

}

int main(int argv, char *argc[]) {

    cerr << "Starting MINA..." << endl;
    Param param(argv,argc);
    cerr << "Getting parameters...";
    ofstream outPa(param.outpath+"param.txt");
    param.getParamInfo(outPa);
    outPa.close();
    cerr << "DONE!" << endl;

    Outcome outcome(param.clinical.c_str());
    vector<string> genenames = readSymbols(param.geneInfo.c_str());
    vector<Profile> profiles;


    for( size_t i = 0 ; i < param.profiles.size() ; ++i ) {
        profiles.push_back( Profile( param.profiles[i].c_str(), 5 ) );
    }
    if( param.method == "network" ) {
        vector<double> thresholds;


        cerr << "Reading profiles and get threshold value..." << endl;
        ofstream outTh(param.outpath+"threshold.txt");

        for( size_t i = 0 ; i < param.profiles.size() ; ++i ) {
            cerr << param.profiles[i] << " ";
            thresholds.push_back( getThreshold( profiles.back(), outcome, param.maxPerm ) );
            cerr << " DONE!" << endl;
            outTh << param.profiles[i] << "\t" << thresholds.back() << endl;
        }
        outTh.close();

        cerr << "Generating networks..." << endl;
        getEdges( genenames, param, profiles, outcome, thresholds );
        cerr <<  " DONE!" << endl;

        cerr << "Finish! You can check the results into output folder!" << endl;
    }
    else if(param.method == "dist") {
        getDist( param, profiles, outcome );
        cerr << "Finish! You can check the results into output folder!" << endl;
    }
    return 0;
}

