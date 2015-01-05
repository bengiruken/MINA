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
#include <random>
#include "Utility.h"
#include "Outcome.h"
#include "Profile.h"
#include "Measure.h"
#include "Param.h"
#include "Network.h"
#include "FrequencyTable.h"

using namespace std;

int main(int argv, char *argc[]) {

	cerr << "Start reading input files." << endl;
	Profile mRNAProfile("C:\\DB\\OV\\mRNA_genomicMatrix_OS", 5);
	Outcome mRNAOutcome("C:\\DB\\OV\\mRNA_outcome_OS.txt");
	auto symbols = readSymbols("C:\\DB\\OV\\mRNA_genomicMatrix_sym");
	cerr << "Reading input files were done." << endl;

	const int nFeatures = mRNAProfile.getNumFeatures();
	FrequencyTable ftAssociation(4);
	FrequencyTable ftInteraction(4);

	showProgress(0, 1, true);
	const long long totalIteration = nFeatures * (nFeatures - 1) / 2;
	long long iteration = 0;
	for (int i = 0; i < (long long)nFeatures; ++i) {
		for (int j = i + 1; j < (long long)nFeatures; ++j) {
			showProgress(++iteration, totalIteration);
			/*
			cerr << symbols[i] << " " << symbols[j];
			*/
			double associationValue = getOutcomeAssociation(
				mRNAProfile[i], mRNAProfile.getNumTypes(),
				mRNAProfile[j], mRNAProfile.getNumTypes(),
				mRNAOutcome);

			double interactionValue = getInteraction(
				mRNAProfile[i], mRNAProfile.getNumTypes(),
				mRNAProfile[j], mRNAProfile.getNumTypes());
			/*
			cerr << " " << associationValue << " " << interactionValue << endl;
			*/
			ftAssociation.put(associationValue);
			ftInteraction.put(interactionValue);
		}
	}
	cerr << endl;

	
	ftAssociation.output("C:\\DB\\OV\\association_debug.txt", 50);
	ftInteraction.output("C:\\DB\\OV\\interaction_debug.txt", 50);
		
	system("PAUSE");
	return 0;
}

