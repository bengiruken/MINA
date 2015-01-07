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
#include <chrono>

#include <ppl.h>
using namespace concurrency;


using namespace std;

int main(int argv, char *argc[]) {

	cerr << "Start reading input files." << endl;
	Profile mRNAProfile("C:\\DB\\OV\\mRNA_genomicMatrix_OS", 5);
	Outcome mRNAOutcome("C:\\DB\\OV\\mRNA_outcome_OS.txt");
	auto symbols = readSymbols("C:\\DB\\OV\\mRNA_genomicMatrix_sym");
	cerr << "Reading input files were done." << endl;

	const int nFeatures = mRNAProfile.getNumFeatures();
	//const int nFeatures = 20;
	FrequencyTable ftAssociation(4);
	FrequencyTable ftInteraction(4);

	auto startTime = chrono::high_resolution_clock::now();

	showProgress(0, 1, true);
	const long long totalIteration = nFeatures * (nFeatures - 1) / 2;
	long long iteration = 0;

	for (int i = 0; i < (long long)nFeatures; ++i) {
		for (int j = i + 1; j < (long long)nFeatures; ++j) {
			showProgress(++iteration, totalIteration);
			double associationValue = getOutcomeAssociation(
				mRNAProfile[i], mRNAProfile.getNumTypes(),
				mRNAProfile[j], mRNAProfile.getNumTypes(),
				mRNAOutcome );

			double interactionValue = getInteraction(
				mRNAProfile[i], mRNAProfile.getNumTypes(),
				mRNAProfile[j], mRNAProfile.getNumTypes());
			ftAssociation.put(associationValue);
			ftInteraction.put(interactionValue);
		}
	}
	auto endTime = std::chrono::high_resolution_clock::now();
	cerr << endl;
	ftAssociation.output("C:\\DB\\OV\\association_real.txt", 30);
	ftInteraction.output("C:\\DB\\OV\\interaction_real.txt", 30);
	cerr << "elpased: " << chrono::duration_cast<chrono::minutes>(endTime - startTime).count() << endl;
	/*
	vector<Outcome> permOutcome(30, mRNAOutcome);
	for (auto p : permOutcome) {
		p.permute();
	}

	random_device rd;
	mt19937 gen(rd());

	double maxAvgAssociation = numeric_limits<double>::min();
	double maxAvgInteraction = numeric_limits<double>::min();

	double associationValue[30];
	double interactionValue[30];

	for (int i = 0; i < (long long)nFeatures; ++i) {

		vector< vector<int> > permProfile(30, mRNAProfile[i]);

		for (auto p : permProfile) {
			shuffle(begin(p), end(p), gen);
		}
		for (int j = i + 1; j < (long long)nFeatures; ++j) {
			showProgress(++iteration, totalIteration);
			for (int k = 0; k < 30; ++k) {

				associationValue[k] = getOutcomeAssociation(
				mRNAProfile[i], mRNAProfile.getNumTypes(),
				mRNAProfile[j], mRNAProfile.getNumTypes(),
				permOutcome[k]);

				interactionValue[k] = getInteraction(
				permProfile[k], mRNAProfile.getNumTypes(),
				mRNAProfile[j], mRNAProfile.getNumTypes());
				
				//cerr << associationValue[k] << " " << interactionValue[k] << endl;
			}

			double avgAssociationValue = accumulate(associationValue, associationValue + 30, 0.0) / 30.0;
			double avgInteractionValue = accumulate(interactionValue, interactionValue + 30, 0.0) / 30.0;
			
			
			ftAssociation.put(avgAssociationValue);
			ftInteraction.put(avgInteractionValue);

			maxAvgAssociation = max(maxAvgAssociation, avgAssociationValue);
			maxAvgInteraction = max(maxAvgInteraction, avgInteractionValue);
		}
	}
	cerr << endl;
	cerr << "random association max = " << maxAvgAssociation << endl;
	cerr << "random interaction max = " << maxAvgInteraction << endl;
	
	ftAssociation.output("C:\\DB\\OV\\association_random.txt", 30);
	ftInteraction.output("C:\\DB\\OV\\interaction_random.txt", 30);
	*/
	system("PAUSE");
	return 0;
}

