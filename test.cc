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

    Network associationNetwork = getAssocationNetwork( profiles, outcome, 0.0, associationThreshold );
    return 0;
}

