#ifndef POLY_FIT_H
#define POLY_FIT_H
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h> // fabs
using namespace std;
class PolyFit
{

public: 
	enum SATUS{
		SUCCESS = 0,
		ERROR = 1
	};

	/*
	* Given a set of 2D points, stored in vectors x and y,
	* finds polynomial coeficients of a polynomial of order
	* pol_order that minimizes the sum of square distances,
	* using LUP-decomposition.
	*
	* INPUT:
	*	- 	Vectors x and y, of equal length corresponding
	*			to a 2D dataset.
	*	- 	Integer pol_order corresponding to the desired order
	*	 		of the polynomial.
	* OUTPUT:
	*	-	Integer, 0 if sucess, 1 if failure.
	*	-	A vector poly_coef, of length pol_order+1,
	* 			with estimated polynomial coeficients. The
	*			coeficients are ordered as:
	*			poly_coef[0] + poly_coef[1] x + ...
	*			 + poly_coef[pol_order] x^pol_order.
	* NOTES:
	*	-	Prints message if matrix is close to degenerate
	*/
	static PolyFit::SATUS polyFit(
		vector<double>&outPolCoef, 
		const vector<double> &inX,
		const vector<double>& inY, 
		const unsigned int& inPolOrder
	);

private:

	/* 
	* INPUT: A - vector of vector, corresponding to a square matrix having dimension N
	*        Tol - small tolerance number to detect failure when the matrix is near degenerate
	* OUTPUT: Matrix A is changed, it contains both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
	*        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1 
	*        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N, 
	*        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S    
	* SOURCE: https://en.wikipedia.org/wiki/LU_decomposition
	*/
	static PolyFit::SATUS LUPDecompose(
		vector<vector<double> >& A,
		const unsigned int& N,
		const double& Tol, 
		vector<unsigned int> &P
	);

	/* Solves the linear system Ax=b, using LUP decomposition.
	*  INPUT: The square matrix A, of size N, and vector b.
	*  OUTPUT: Vector x, the solution of Ax = b.
	*  SOURCE: https://en.wikipedia.org/wiki/LU_decomposition
	*/
	static PolyFit::SATUS LUPSolve(
		vector<double>& x,
		vector<vector<double> >& A,
		const vector<double>& b,
		const unsigned int& N
	);
};


#endif
