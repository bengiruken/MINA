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

    const int maxPerm = 30;
    Outcome outcome("/home/hhjeong/ws/OV/clinical.txt");
    vector<string> genenames = readSymbols("/home/hhjeong/ws/OV/sym.txt");

    Profile profiles( "/home/hhjeong/ws/OV/mRNA.txt", 5 );
    double interactionThreshold = getInteractionNetworkThreshold( profiles, maxPerm );
    double associationThreshold = getOutcomeAssociationNetworkThreshold( profiles, outcome, maxPerm );

    cerr << interactionThreshold << endl;
    cerr << associationThreshold << endl;


    double alpha[] = { 0.0, 0.2, 0.4, 0.6, 0.8 };

    for( int i = 0 ; i < 5 ; ++i ) {
        cerr << "[DEBUG] alpha = " << alpha[i] << endl;
        Network associationNetwork = getAssocationNetwork( profiles, outcome, 
            associationThreshold, alpha[i] );
        Network interactionNetwork = getInteractionNetwork( profiles, 
            interactionThreshold, alpha[i] );
        
        char fileName[128];

        sprintf(fileName,"network/mRNA_association_network_%.1f.txt", alpha[i]);
        saveNetwork( associationNetwork, fileName );
        sprintf(fileName,"network/mRNA_interaction_network_%.1f.txt", alpha[i]);
        saveNetwork( interactionNetwork, fileName );


        Network pureAssociationNetwork = getDifferenceNetwork( associationNetwork, interactionNetwork );
        sprintf(fileName,"network/mRNA_pure_network_%.1f.txt", alpha[i]);
        saveNetwork( pureAssociationNetwork, fileName );
    }
    return 0;
}

