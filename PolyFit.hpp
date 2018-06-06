#ifndef POLY_FIT_H
#define POLY_FIT_H
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h> // fabs
using namespace std;


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
int polyFit(vector<double>&outPolCoef, 
    const vector<double> &inX,
	const vector<double>& inY, 
	const unsigned int& inPolOrder);
#endif
