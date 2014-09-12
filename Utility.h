#ifndef __FILEHANDLER__
#define __FILEHANDLER__
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;

template <typename T>
vector< vector<T> > read(  const char *filename ) {
	ifstream inp(filename);
	string line;	

	vector< vector<T> > ret;
	int i = 0;
	while( getline(inp,line) ) {
		istringstream sin(line);

		string genename;
		vector<T> row;
		T value;
		sin >> genename;
		while( sin >> value ) {
			row.push_back(value);
		}

		ret.push_back(row);
	}
	inp.close();	
	return ret;
}

vector< string > readSymbols( const char *filename ) {
	vector< string > ret;
	ifstream inp(filename);
	string tok;
	while( inp >> tok ) {
		ret.push_back(tok);
	}
	return ret;
}

template <typename T>
double getMedian( vector<T> seq ) {
	sort( seq.begin(), seq.end() );
	const size_t size = seq.size();
	if( size % 2 == 1 ) {
		return seq[size/2];
	}
	else {
		return 0.5 * ( seq[size/2-1] + seq[size/2] );
	}
}

#endif
