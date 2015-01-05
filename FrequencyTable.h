#ifndef __FREQUENCY_TABLE__
#define __FREQUENCY_TABLE__

#include <vector>
#include <utility>
#include <algorithm>
using namespace std;

// for [0,1] inteval
class FrequencyTable {
private:
	int order;
	vector<long long> freq;
public:
	FrequencyTable(const int precision = 4) {
		// if the precision is set 4 then range which considered in this class is [0.0001,1.0000]
		order = 1;
		for (int i = 0; i < precision; ++i) order *= 10;
		freq = vector<long long>(order + 1, 0);

	}

	void clear() {
		fill(freq.begin(), freq.end(), 0);
	}

	void put(const double &x) {
		freq[int(x*order)]++;
	}

	int getMinPos() {
		for (int i = 0; i < (int)freq.size(); ++i) {
			if (freq[i] >0) {
				return i;
			}
		}
		return -1;
	}

	int getMaxPos() {
		for (int i = (int)freq.size() - 1; i >= 0;  --i) {
			if (freq[i] > 0) {
				return i;
			}
		}
		return -1;
	}

	vector< pair<double, long long> > freq2bin(int bin) {
		int mini = getMinPos();
		int maxi = getMaxPos();

		const int width = (int)ceil((maxi - mini) * 1.0 / bin);
		vector< pair<double,long long> > ret;
		for (int i = mini; i <= maxi; i += width) {
			long long cumsum = 0;
			for (int j = i; j < min(i + width, maxi+1); ++j) {
				cumsum += freq[j];
			}
			ret.push_back(make_pair((double)i / order, cumsum));
		}

		return ret;
	}
};

#endif