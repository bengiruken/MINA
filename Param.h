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
    string method;
    string outpath;
	string measure;

    int maxPerm;
    double alpha;
    double distLo, distHi;
	int nbins;

    Param(int argv, char *argc[]) {
        distLo = 0.0;
        distHi = 1.0;
		nbins = 30;
        maxPerm = 30;
        outpath = "output\\";
		measure = "omi";
        for( int i = 1 ; i < argv ; ++i ) {
            if( strcmp( argc[i], "-s") == 0 ) {
                geneInfo = argc[++i];
            }
            else if( strcmp( argc[i], "-ip") == 0 ) {
                profiles.push_back(argc[++i]);
            }
            else if( strcmp( argc[i], "-io") == 0 ) {
                clinical = argc[++i];
            }
            else if( strcmp( argc[i], "-o") == 0 ) {
                outpath = argc[++i];
            }
            else if( strcmp(argc[i], "-perm") == 0 ) {
                sscanf(argc[++i], "%d", &maxPerm);
            }
            else if( strcmp(argc[i], "-alpha") == 0 ) {
                sscanf(argc[++i], "%lf", &alpha);
            }
            else if( strcmp(argc[i], "-dlo") == 0 ) {
                sscanf(argc[++i], "%lf", &distLo);
            }
            else if( strcmp(argc[i], "-dhi") == 0 ) {
                sscanf(argc[++i], "%lf", &distHi);
            }
			else if( strcmp(argc[i], "-dbin") == 0 ) {
				sscanf(argc[++i], "%d", &nbins);
			}
            else if( strcmp(argc[i],"dist") == 0 || strcmp(argc[i],"network") == 0 ) {
                method = argc[i];
            }
			else if( strcmp(argc[i], "-m") == 0 ) {
				measure = argc[++i];
			}
            else {
                errorMsg( string(argc[i]) + " is not a valid parameter" );
            }
        }
        
    }

    void getParamInfo( ostream &out ) {
            out << "Informations" << endl;
            out << "Start time : " << currentDateTime() << endl;
            out << "method : " << method << endl;
            out << "geneInfo : " << geneInfo << endl;
            out << "profiles : " << join( profiles, "," ) << endl;
            out << "clinical : " << clinical << endl;
            out << "maxPerm : " << maxPerm << endl;
            out << "alpha :" << alpha << endl;
            out << "outpath : " << outpath << endl;
    }
};

#endif
