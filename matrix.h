#ifndef Matrix_h
#define Matrix_h

#include <vector>

using namespace std;

class Eigenpair
{
	public:
		Eigenpair(double val, double * vec);
		double value;
		double * vector;
};

class Matrix 
{
	private:
	 	double ** mat;
	 	int rows, cols;
	 	Matrix();
	 	void getCharPol(int rowNum, vector<int> exCols, double * pol);
	 	void getEigenvector(double val, double * iVec);
	public:
		Matrix(int rs, int cls);
		void set(int row, int col, double val);
		double get(int row, int col);
		void add(int row, int col, double val);
		void print();
		void rowReduce();
		vector<Eigenpair> getEigenpairs();
		Matrix copy();
};

#endif
