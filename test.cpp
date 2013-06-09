#include <iostream>
#include <cmath>
#include "polynomial.h"
#include "matrix.h"

#define NUM_POINTS 10

using namespace std;

class Point
{
	public:
		Point (double, double);
		double x, y;
};

Point::Point(double _x, double _y)
{
	x = _x;
	y = _y;
}

void printM(double ** mat, int rows, int cols)
{
	for (int j = 0; j < rows; j++)
	{
		for (int i = 0; i < cols; i++)
		{
			cout << *(*(mat + j) + i) << " ";
		}
		cout << endl;
	}
}

void solveQuad(double a, double b, double c, double * outp)
{
	*(outp++) = (-b + sqrt(b * b - 4 * a * c)) / 2 * a;
	*outp = (-b - sqrt(b * b - 4 * a * c)) / 2 * a;
}


int main() 
{	
	/*double prepMat[3][3] = { { 4, 6, 0}, {4, 0, 1}, {2, 3, 0} };
	
	double * mat[3] = {prepMat[0], prepMat[1], prepMat[2]};
	
	printM(mat, 3, 3);
	cout << "-----" << endl;
	
	rowReduce(mat, 3, 3);
	
	printM(mat, 3, 3);
	
	double nullVec[3];
	findNullVector(mat, 3, 3, nullVec);
	double * vecSet[1] = {nullVec};
	cout << "eigen values: ";
	printM(vecSet, 1, 3);*/
	/*
	
	Point ps [NUM_POINTS] = {Point(-1, -1), Point(1, 1), Point(-2, -2), Point(2, 2), 
		Point(-4, -3), Point(4, 3), Point(-4, -5), Point(4, 5), Point(-5, -5), Point(5, 5)};
	double mat[2][2];
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			mat[i][j] = 0;
		}
	} 
	for (int i = 0; i < NUM_POINTS; i++)
	{
		mat[0][0] += ps[i].x * ps[i].x;
		mat[0][1] += ps[i].x * ps[i].y;
		mat[1][0] += ps[i].y * ps[i].x;
		mat[1][1] += ps[i].y * ps[i].y;
	}
	
	cout << "matrix: " << endl;
	double * ttt[2] = { mat[0], mat[1] };
	printM(ttt, 2, 2);
	
	double vals[2];
	solveQuad(1, -mat[0][0] - mat[1][1], mat[0][0] * mat[1][1] - 
		mat[0][1] * mat[1][0], vals);
	//cout << "eigenvalues: " << endl;
	//cout << vals[0] << ", " << vals[1] << endl;
	//cout << "eigen vectors: " << endl;
	for (int i = 0; i < 2; i++)
	{
		cout << "eigen value " << vals[i] << " vector: ";
		double eigVec[2];
		getEigenvector(mat, vals[i], eigVec);
		cout << eigVec[0] << ", " << eigVec[1] << endl;
	}
	cout << "hello" << endl;
	
	double polCfs[3] = { -1, 0, 1 };
	Polynomial pol(polCfs, 3);
	vector<double> roots = pol.getRoots();
	for (vector<double>::iterator it = roots.begin(); it != roots.end(); ++it)
    cout << *it << endl;*/
    
  Matrix mat(3, 3);
  mat.set(2, 2, 2.132);
  mat.set(2, 1, 21.32);
  mat.set(2, 0, 213.2);
  mat.set(1, 1, 2132.13);
  mat.set(1, 2, 21.2);
  //cout << mat.get(0, 0) << " " << mat.get(2, 2) << endl;
  mat.print();
}
