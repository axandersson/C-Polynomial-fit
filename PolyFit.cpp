#include "PolyFit.hpp"


PolyFit::SATUS PolyFit::LUPDecompose(
	vector<vector<double> >& A,
	const unsigned int& N,
	const double& Tol,
	vector<unsigned int> &P
)
{
	unsigned int i, j, k, imax;
	double maxA, absA;

	for (i = 0; i <= N; i++)
		P[i] = i; 

	for (i = 0; i < N; i++) 
	{
		maxA = 0.0;
		imax = i;

		for (k = i; k < N; k++)
			if ((absA = fabs(A[k][i])) > maxA) 
			{
				maxA = absA;
				imax = k;
			}

		if (maxA < Tol) return SATUS::ERROR; //failure, matrix is degenerate

		if (imax != i) 
		{
			j = P[i];
			P[i] = P[imax];
			P[imax] = j;

			A[i].swap(A[imax]);
			P[N]++;
		}

		for (j = i + 1; j < N; j++) 
		{
			A[j][i] /= A[i][i];

			for (k = i + 1; k < N; k++)
				A[j][k] -= A[j][i] * A[i][k];
		}
	}
	return SATUS::SUCCESS;  //decomposition done 
}


PolyFit::SATUS PolyFit::LUPSolve(
	vector<double>& x,
	vector<vector<double> >& A,
	const vector<double>& b,
	const unsigned int& N
)
{
	const double tol = 1e-12;
	vector<unsigned int> P(N + 1,0);

	if(LUPDecompose(A, N, tol, P) == SATUS::SUCCESS)
	{
		for (unsigned int i = 0; i < N; i++) 
		{
			x[i] = b[P[i]];

			for (unsigned int k = 0; k < i; k++)
				x[i] -= A[i][k] * x[k];
		}

		for (int i = N - 1; i >= 0; i--) 
		{
			for (unsigned int k = i + 1; k < N; k++)
				x[i] -= A[i][k] * x[k];

			x[i] = x[i] / A[i][i];
		}
	}
	else
	{

		return PolyFit::SATUS::ERROR;
	}
	return PolyFit::SATUS::SUCCESS;
}

PolyFit::SATUS PolyFit::polyFit(
	vector<double>& outPolCoef, 
    const vector<double>& inX,
	const vector<double>& inY, 
	const unsigned int& inPolOrder
)
{
	const unsigned int nX = inX.size();
	const unsigned int nCoef = inPolOrder + 1;
	const unsigned int nPower = nCoef + inPolOrder;

	outPolCoef.resize(nCoef,0.0);

	// The normal matrices.
	// (X' X) * poly_coef = (X' y).
	vector<vector<double> > XX(nCoef, vector<double>(nCoef, 0));
	vector<double > XY(nCoef, 0);
	
	vector<double > xxSum(nPower, 0);

	// Calculating all powers necessary for constructing
	// the normal matrix X' X, and X' y.
	for (unsigned int i = 0; i < nX; i++)
	{
		double accumulatedPower = 1;
		for (unsigned int j = 0; j < nPower; j++)
		{
			if (j < nCoef)
			{
				XY[j] += inY[i] * accumulatedPower;
			}
			xxSum[j] += accumulatedPower;
			accumulatedPower *= inX[i];
		}
	}
	for (unsigned int i = 0; i < nCoef; i++)
	{
		for (unsigned int j = 0; j < nCoef; j++)
		{
			XX[i][j] = xxSum[i + j];
			XX[j][i] = XX[i][j];
		}
	}

	const SATUS errCode = LUPSolve(outPolCoef,XX, XY, nCoef);
	return errCode;
}
