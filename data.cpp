#include "data.h"

using namespace std;

vector<Eigenpair> Data::getPairs(double * points, int numPts, int dim)
{
	//get avg for each coord
	double avgVec[dim];
	for (int i = 0; i < dim; i++)
		avgVec[i] = 0;
	for (int ptInd = 0; ptInd < numPts; ptInd++)
	{
		for (int dmInd = 0; dmInd < dim; dmInd++)
		{
			avgVec[dmInd] += *(points + ptInd * dim + dmInd) / numPts;
		}
	}
		
	//compute covariance matrix
	Matrix covM(dim, dim);
	for (int ptInd = 0; ptInd < numPts; ptInd++)
	{
		for (int i = 0; i < dim; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				covM.add(i, j, (*(points + ptInd * dim + i) - avgVec[i]) *
					(*(points + ptInd * dim + j) - avgVec[j]));
			}
		}
	}
	
	//return eigen pairs for covariance matrix
	return covM.getEigenpairs();
}
