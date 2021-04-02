#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* uses Householder QR to compute LS solution (see Golub & Loan) */
double *hhls(double **A, double *b, int m, int n) {

	int i, j, p, q; /* counters */
	double *v, s, sigma, *beta, mu; /* Householder vector, a scalar, sigma, Householder beta and mu */
	
	/* allocate memory for the Householder vector */
	if ((v = (double *) malloc(m*sizeof(double))) == NULL)
		return NULL;

	/* allocate memory for the Householder beta */
	if ((beta = (double *) malloc(n*sizeof(double))) == NULL)
		return NULL;

	/* rewrite A with its Householder QR decomposition */
	for (j = 0; j < n; j++) {

		/* compute householder vector and beta */
		sigma = 0.;
		for (i = j + 1; i < m; i++)
			sigma += A[i][j] * A[i][j];
		v[j] = 1.;
		for (i = j + 1; i < m; i++)
			v[i] = A[i][j];
		if (sigma == 0.)
			beta[j] = 0.;
		else {
			mu = sqrt(A[j][j] * A[j][j] + sigma);
			if (A[j][j] <= 0.)
				v[j] = A[j][j] - mu;
			else
				v[j] = -sigma / (A[j][j] + mu);
			beta[j] = 2 * v[j] * v[j] / (sigma + v[j] * v[j]);
			for (i = m - 1; i >= j; i--)
				v[i] = v[i] / v[j];
		}

		/* apply the Householder transformation to A */
		for (q = j; q < n; q++) {

			s = 0.;
			for (p = j; p < m; p++)
				s += v[p] * A[p][q];

			for (p = j; p < m; p++)
				A[p][q] -= beta[j] * v[p] * s;
		}

		/* copy the Householder vector into A */
		if (j < m - 1)
			for (i = j + 1; i < m; i++)
				A[i][j] = v[i];
	}

	/* apply Householder transformations to b */
	for (j = 0; j < n; j++) {

		v[j] = 1.;
		for (i = j + 1; i < m; i++)
			v[i] = A[i][j];

		s = 0.;
		for (i = j; i < m; i++)
			s += v[i] * b[i];

		for (i = j; i < m; i++)
			b[i] -= beta[j] * v[i] * s;

	}

	/* backsolve */
	b[n - 1] /= A[n - 1][n - 1];
	for (i = n - 2; i >= 0; i--) {
		b[i] /= A[i][i];
		for (j = i + 1; j < n; j++)
			b[i] -= A[i][j] * b[j] / A[i][i];
	}

	free(v); 
	free(beta);

	return b;
}

