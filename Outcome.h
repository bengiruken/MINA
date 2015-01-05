#ifndef __OUTCOME__
#define __OUTCOME__
#include <vector>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <random>

using namespace std;

class Outcome {
private:
	random_device rd;
	mt19937 gen;
	size_t numTypes;
	vector<int> status;	
	vector<int> numSubjects;

public:
	int readFromFile( const char *fileName ) {
		ifstream inp(fileName);
		
		int num;
		status.clear();
		while( inp >> num ) {
			status.push_back(num);
		}
		int mini = *min_element(status.begin(),status.end());
		int maxi = *max_element(status.begin(),status.end());
		assert( 0 == mini );

		inp.close();

		return maxi-mini+1;
	}

	Outcome( const char *fileName ) {
		numTypes = readFromFile( fileName );
		numSubjects.resize( numTypes, 0 );
		for( vector<int>::iterator it = status.begin() ; it != status.end() ; ++it ) {
			numSubjects[*it]++;
		}
		gen = mt19937(rd());
	}

	Outcome( const Outcome &other ) {
		numTypes = other.numTypes;
		numSubjects = other.numSubjects;
		status = other.status; 
	}

	size_t getNumTypes() const {
		return numTypes;	
	}

	int getNumSubjects( int type ) const {
		return numSubjects[type];
	}

	size_t size() const {
		return status.size();
	}
	
	int operator[]( const int &index )  {
		return status[index];
	}

	void permute() {
		shuffle( status.begin(), status.end(), gen );
	}
};

#endif
