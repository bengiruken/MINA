#ifndef __PROFILE__
#define __PROFILE__
#include "Utility.h"
#include <vector>
#include <cassert>
#include <limits>
#include <map>

using namespace std;

class Profile { 

private:
    int numTypes;
    vector< vector<int> > dat;

    vector< int > binning( const vector<int> &row ) {
        vector< int > ret;

        int mini = *min_element( row.begin(), row.end() );
        for( size_t i = 0 ; i < row.size() ; ++i ) {
            ret.push_back( row[i] - mini );
        }

		numTypes = *max_element(ret.begin(),ret.end()) - mini + 1;
        return ret;
    }

	vector< int > binning(const vector<string> &row) {
		vector< int > ret;

		int number = 0;
		map< string, int > str2int;

		for ( auto x : row) {
			if (str2int.count(x) == 0) {
				str2int[x] = number++;
			}
		}

		for ( auto x : row) {
			ret.push_back(str2int[x]);
		}

		numTypes = number;

		return ret;
	}

    vector< int > binning( const vector<double> &row, int size = 5 ) {
        vector< int > ret;
        auto mini = *min_element( row.begin(), row.end() );
        auto maxi = *max_element( row.begin(), row.end() );

        double scale = maxi - mini;
        for( auto x : row ) {
            double val = x - mini;
            int normVal = min( size-1, (int)(val / scale * size) );
            ret.push_back(normVal);
        }

        return ret;
    }

public:
    Profile( const char *filename, int binSize = 5 ) {
        vector< vector<double> > src = read<double>(filename);
        int mini = numeric_limits<int>::max();
		int maxi = numeric_limits<int>::min();

        for( size_t i = 0 ; i < src.size() ; ++i ) {
            dat.push_back( binning( src[i], binSize ) );
            for( size_t j = 0 ; j < dat.back().size() ; ++j ) {
                mini = min( mini, dat[i][j] );
                maxi = max( maxi, dat[i][j] );
            }	
        }
        numTypes = maxi - mini + 1;
        assert( 0 == mini );
    }

    size_t getNumTypes() {
        return numTypes;
    }

    size_t getNumFeatures() {
        return dat.size();
    }
    size_t getNumRows() {
        return dat.front().size();
    }

    vector<int>& operator[] ( int index ) {
        return dat[index];
    }

};
#endif
