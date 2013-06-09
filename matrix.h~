#ifndef Matrix_h
#define Matrix_h

#include <vector>

class Eigenpair
{
	public:
		double value;
		double * vector;
};

class Matrix 
{
	private:
	 	double ** mat;
	 	int rows, cols;
	 	Matrix();
	 	void getCharPol(int rowNum, vector<int> cols, double * pol);
	public:
		Matrix(int rs, int cls);
		void set(int row, int col, double val);
		double get(int row, int col);
		void print();
		void rowReduce();
		vector<Eigenpair> getEigenpairs();
		Matrix copy();
};

#endif
