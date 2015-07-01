#ifndef __FREQUENCY_TABLE__
#define __FREQUENCY_TABLE__

#include <vector>
#include <utility>
#include <algorithm>
using namespace std;

class FrequencyTable {
private:
	double lo, hi;
	double interval;
	int bins;
	vector<long long> freq;
public:
	FrequencyTable(const double _lo = 0, const double _hi = 1, const int _bins = 30) {
		lo = _lo;
		hi = _hi;
		bins = _bins;
		interval =  (hi-lo)/bins;
		freq = vector<long long>(bins, 0);	
	}

	void clear() {
		fill(freq.begin(), freq.end(), 0);
	}

	void put(const double &x) {
		freq[ min(bins-1, int(x/interval)) ]++;
	}

	void output(const char *outName ) {
		ofstream oup(outName);
		for( int i = 0 ; i < bins ; ++i ) {
			oup << lo + i*interval << " " << freq[i] << endl;
		}
		oup.close();
	}

};

#endif
