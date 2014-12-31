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

    const double alpha[] = { 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2 };

    const double associationThreshold = 0.07;
    const double interactionThreshold = 0.04;

    const Network interactionNetwork = 
        loadNetwork( "network/mRNA_interaction_network_0.0.txt" );

    const Network associationNetwork =
        loadNetwork( "network/mRNA_association_network_0.0.txt" );

    cerr << "information of interactionNetwork" << endl;
    getNetworkInformation( interactionNetwork );
    cerr << "information of associationNetwork" << endl;
    getNetworkInformation( associationNetwork );


    for( int i = 1 ; i < 7 ; ++i ) {
        cerr << "[DEBUG] alpha = " << alpha[i] << endl;
        Network filteredAssociationNetwork = filterNetwork( associationNetwork, 
            (1+alpha[i]) * associationThreshold );
        Network filteredInteractionNetwork = filterNetwork( interactionNetwork,
            (1+alpha[i]) * interactionThreshold );

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

