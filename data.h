#ifndef data_h
#define data_h

#include "matrix.h"

#include <vector>

using namespace std;

class Data
{
	public:
		vector<Eigenpair> getPairs(double * points, int numPts, int dim);
};

#endif
