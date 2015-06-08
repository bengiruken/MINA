#ifndef __PARAM__
#define __PARAM__

#include <cstdio>
#include <fstream>
#include <string>
#include "Utility.h"
using namespace std;

struct Param {
    string geneInfo;
    vector<string> profiles;
    string clinical;
    string separation;
    int maxPerm;
    vector<double> alpha;

    Param(int argv, char *argc[]) {
        string line; 
        getline( inp, line );
        while ( getline( inp, line ) ) {
            vector<string> toks = tokenization( line );

            if( toks.front() == "geneinfo:" ) {
                geneInfo = toks.back(); 
            }
            else if( toks.front() == "profiles:" ) {
                for( size_t i = 1 ; i < toks.size() ; ++i ) {
                    profiles.push_back(toks[i]);
                }
            }
            else if( toks.front() == "clinical:" ) { 
                clinical = toks.back();
            }
            else if( toks.front() == "separation:" ) {
                separation = toks.back();
            }
            else if( toks.front() == "maxPerm:" ) {
                sscanf(toks.back().c_str(), "%d", &maxPerm );
            }
            else if( toks.front() == "alpha:" ) {
                for( size_t i = 1 ; i < toks.size() ; ++i ) {
                    double val;
                    sscanf(toks[i].c_str(),"%lf", &val);
                    alpha.push_back(val);
                }
            }
            else {
                errorMsg( toks.front() + " is not a valid parameter" );
            }
        }
        
        inp.close();
    }

    void getParamInfo( ostream &out ) {
            out << "Informations" << endl;
            out << "Start time : " << currentDateTime() << endl;
            out << "geneInfo : " << geneInfo << endl;
            out << "profiles : " << join( profiles, "," ) << endl;
            out << "clinical : " << clinical << endl;
            out << "separation : " << separation << endl;
            out << "maxPerm : " << maxPerm << endl;
            out << "alpha :";
            for( size_t i = 0 ; i < alpha.size() ; ++i ) {
                out << " " << alpha[i];
            }
            out << endl;
    }
};

#endif
