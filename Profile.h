#ifndef __PROFILE__
#define __PROFILE__
#include "Utility.h"
#include <vector>
#include <cassert>
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

        return ret;
    }

    vector< int > binning( const vector<double> &row, int size = 5 ) {
        vector< int > ret;
        double mini = *min_element( row.begin(), row.end() );
        double maxi = *max_element( row.begin(), row.end() );

        double scale = maxi - mini;
        for( size_t i = 0 ; i < row.size() ; ++i ) {
            double val = row[i] - mini;
            int normVal = min( size-1, (int)(val / scale * size) );
            ret.push_back(normVal);
        }

        return ret;
    }

public:
    Profile( const char *filename, int binSize = 5 ) {
        vector< vector<double> > src = read<double>(filename);
        int mini = 2147483647;
        int maxi = -2147483648;

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
