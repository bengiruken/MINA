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
#include "Network.h"

using namespace std;

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
        cerr << param.profiles[i] << " ";
        thresholds.push_back( getOutcomeAssociationNetworkThreshold( profiles.back(), outcome, param.maxPerm ) );
        cerr << " DONE!" << endl;
        outTh << param.profiles[i] << "\t" << thresholds.back() << endl;
    }
    outTh.close();

    cerr << "Generating networks..." << endl;
    // getAssocationNetwork( genenames, param, profiles, outcome, thresholds );
    cerr <<  " DONE!" << endl;

    cerr << "Finish! You can check the results into output folder!" << endl;
    return 0;
}

