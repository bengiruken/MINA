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

    Profile profile( "/home/hhjeong/ws/OV/mRNA.txt", 5 );
    double interactionThreshold = getInteractionNetworkThreshold( profile, maxPerm );
    double associationThreshold = getOutcomeAssociationNetworkThreshold( profile, outcome, maxPerm );

    cerr << interactionThreshold << endl;
    cerr << associationThreshold << endl;

    double alpha[] = { 0.0, 0.2, 0.4, 0.6, 0.8 };

    Network associationNetwork = getAssocationNetwork( profile, outcome, 
            associationThreshold );
    Network interactionNetwork = getInteractionNetwork( profile, 
            interactionThreshold );

    for( int i = 0 ; i < 5 ; ++i ) {
        cerr << "[DEBUG] alpha = " << alpha[i] << endl;
        Network filteredAssociationNetwork = filterNetwork( associationNetwork, 
            alpha[i] * associationThreshold );
        Network filteredInteractionNetwork = filterNetwork( interactionNetwork,
            alpha[i] * interactionThreshold );

        char fileName[128];

        sprintf(fileName,"network/mRNA_association_network_%.1f.txt", alpha[i]);
        saveNetwork( filteredAssociationNetwork, fileName );
        sprintf(fileName,"network/mRNA_interaction_network_%.1f.txt", alpha[i]);
        saveNetwork( filteredInteractionNetwork, fileName );

        Network pureAssociationNetwork = getDifferenceNetwork( 
                filteredAssociationNetwork, filteredInteractionNetwork );

        sprintf(fileName,"network/mRNA_pure_network_%.1f.txt", alpha[i]);
        saveNetwork( pureAssociationNetwork, fileName );
    }
    return 0;
}

