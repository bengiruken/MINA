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

	const string path = "D:\\DB\\OV\\";
	const string profileName[] = { "mRNA", "CNA", "METH" };
	const string symbolFileName = "D:\\DB\\OV\\sym.txt";
	const string outcomeFileName = "D:\\DB\\OV\\clinical.txt";

	Outcome outcome(outcomeFileName.c_str());
	auto symbols = readSymbols(symbolFileName.c_str());

	for (auto p : profileName) {
		string profileFileName = path + p + ".txt";
		cout << "current: " << p << endl;
		auto startTime = chrono::high_resolution_clock::now();
		Profile pf(profileFileName.c_str(), 5);
		auto endTime = std::chrono::high_resolution_clock::now();
		cout << "elpased: " << chrono::duration_cast<chrono::seconds>(endTime - startTime).count() << "sec" << endl;

		cout << "start calculation" << endl;

		const int nFeatures = pf.getNumFeatures();
		showProgress(0, 1, true);
		FrequencyTable ftAssociation(4);
		FrequencyTable ftInteraction(4);

		vector<Outcome> permOutcome(30, outcome);
		for (auto p : permOutcome) {
			p.permute();
		}

		random_device rd;
		mt19937 gen(rd());

		double maxAvgAssociation = numeric_limits<double>::min();
		double maxAvgInteraction = numeric_limits<double>::min();

		double associationValue[30];
		double interactionValue[30];

		long long iteration = 0;
		long long totalIteration = nFeatures * (nFeatures - 1) / 2;
		startTime = chrono::high_resolution_clock::now();
		for (int i = 0; i < (long long)nFeatures; ++i) {

			vector< vector<int> > permProfile(30, pf[i]);

			for (auto p : permProfile) {
				shuffle(begin(p), end(p), gen);
			}
			for (int j = i + 1; j < (long long)nFeatures; ++j) {
				showProgress(++iteration, totalIteration);

#pragma loop(hint_parallel(0))
				for (int k = 0; k < 30; ++k) {

					associationValue[k] = getOutcomeAssociation(
						pf[i], pf.getNumTypes(),
						pf[j], pf.getNumTypes(),
						permOutcome[k]);

					interactionValue[k] = getInteraction(
						permProfile[k], pf.getNumTypes(),
						pf[j], pf.getNumTypes());

				}

				double avgAssociationValue = accumulate(associationValue, associationValue + 30, 0.0) / 30.0;
				double avgInteractionValue = accumulate(interactionValue, interactionValue + 30, 0.0) / 30.0;


				ftAssociation.put(avgAssociationValue);
				ftInteraction.put(avgInteractionValue);

				maxAvgAssociation = max(maxAvgAssociation, avgAssociationValue);
				maxAvgInteraction = max(maxAvgInteraction, avgInteractionValue);
			}
		}
		endTime = std::chrono::high_resolution_clock::now();
		cout << "elpased: " << chrono::duration_cast<chrono::seconds>(endTime - startTime).count() << "sec" << endl;
		cout << endl;
		cout << "random association max = " << maxAvgAssociation << endl;
		cout << "random interaction max = " << maxAvgInteraction << endl;

		ftAssociation.output((path+"association_"+p+".txt").c_str(), 30);
		ftInteraction.output((path + "interaction_" + p + ".txt").c_str(), 30);
	}

	system("PAUSE");
	return 0;
}

