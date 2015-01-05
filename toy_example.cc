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
	FrequencyTable ft(4);

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dist;

	for (int i = 0; i < (int)1e5; ++i) {
		ft.put(dist(gen));
	}

	cerr << "minPos: " << ft.getMinPos() << endl;
	cerr << "maxPos: " << ft.getMaxPos() << endl;

	auto ret = ft.freq2bin(30);

	for (auto x : ret) {
		cerr << x.first << "~: " << x.second << endl;
	}

	auto vec = [&ret](int i){ return ret[i].second; };
	cerr << accumulate(ret.begin(), ret.end(), 0,
		[](long long sum, pair<double, long long> x){ 
		return sum + x.second; 
	}) << endl;
	system("PAUSE");
	return 0;
}

