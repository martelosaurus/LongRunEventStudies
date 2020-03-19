#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ecf.h" 
#include "ret.h" 
#include "tio.h"
#include "tls.h" 
#include "tps.h" 
#include "ttc.h"
#include "tuf.h"

/* R^2 Suite */
extern double **A_ls; /* regressor matrix for least-squares */
	extern double **A_temp; /* temporary 'A_ls' */
extern double *b_ls; /* regressand matrix for least-squares */
	extern double *b_temp; /* temporary 'b_ls' */
extern double *v_hh; /* Householder vector */
extern double *b_hh; /* Householder beta */
extern double *b_hat; /* fitted right-hand-side */

/* uses Householder QR to compute LS solution (see Golub & Loan) */
double *hhls(void) {

	short i, j, p, q; /* counters */
	double *v, s, sigma, *beta, mu; /* Householder vector, a scalar, sigma, Householder beta and mu */

	/* copy A_temp, b_temp into A_ls, b_ls, resp. */
	for (i = 0; i < ESTW; i++) {
		for (j = 0; j < MNMC; j++)
			A_ls[i][j] = A_temp[i][j];
		b_ls[i] = b_temp[i];
	}
	
	/* [BEGIN] rewrite A_ls with its Householder QR decomposition */
	for (j = 0; j < MNMC; j++) {
	
		/* [BEGIN] compute householder vector and beta */
		
		sigma = 0.;
		for (i = j+1; i < ESTW; i++)
			sigma += A_ls[i][j]*A_ls[i][j];
		
		v_hh[j] = 1.;
		for (i = j+1; i < ESTW; i++)
			v_hh[i] = A_ls[i][j];

		if (sigma == 0.)
			b_hh[j] = 0.;
		else {
			mu = sqrt(A_ls[j][j]*A_ls[j][j]+sigma);
			if (A_ls[j][j] <= 0.)
				v_hh[j] = A_ls[j][j]-mu;
			else
				v_hh[j] = -sigma/(A_ls[j][j]+mu);
			b_hh[j] = 2*v_hh[j]*v_hh[j]/(sigma+v_hh[j]*v_hh[j]);
			for (i = ESTW-1; i >= j; i--)
				v_hh[i] = v_hh[i]/v_hh[j];
		}
		
		/* [END] */
			
		for (q = j; q < MNMC; q++) {
			
			s = 0.;
			for (p = j; p < ESTW; p++)
				s += v_hh[p]*A_ls[p][q];
				
			for (p = j; p < ESTW; p++)
				A_ls[p][q] -= b_hh[j]*v_hh[p]*s;
		}

		if (j < ESTW-1) 
			for (i = j+1; i < ESTW; i++)
				A_ls[i][j] = v_hh[i];
	}
	
	/* [END] */
	
	/* [BEGIN] apply Householder transformations to right-hand-side */
	for (j = 0; j < MNMC; j++) {
		
		v_hh[j] = 1.;
		for (i = j+1; i < ESTW; i++)
			v_hh[i] = A_ls[i][j];
			
		s = 0.;
		for (i = j; i < ESTW; i++)
			s += v_hh[i]*b_ls[i];
			
		for (i = j; i < ESTW; i++)
			b_ls[i] -= b_hh[j]*v_hh[i]*s;
	
	}
	
	/* [END] */
	
	/* [BEGIN] backsolve */
	b_ls[MNMC-1] /= A_ls[MNMC-1][MNMC-1];
	for (i = MNMC-2; i > -1; i--) {
		b_ls[i] /= A_ls[i][i];
		for (j = i+1; j < MNMC; j++)
			b_ls[i] -= A_ls[i][j]*b_ls[j]/A_ls[i][i];
	}
		
	/* [END] */
	
	return b_ls;
}

/* compute R^2 */
double calc_rsqr(long evnt) {

	short i, j; /* counters */
	double *x_ls, RSS; /* LS solution and RSS */
	static double b_bar, TSS; /* sample mean and TSS */
	static long prev_evnt = -1; /* previous event */
	
	if (evnt != prev_evnt) { /* if the current event is different from the previous event... */
	
		b_bar = TSS = 0.; /* ...initialize... */
	
		/* ...compute the sample mean... */
		for (i = 0; i < ESTW; i++) 
			b_bar += b_temp[i]/ESTW;
		
		/* ...compute the TSS... */
		for (i = 0; i < ESTW; i++)
			TSS += (b_temp[i]-b_bar)*(b_temp[i]-b_bar);
			
		prev_evnt = evnt; /* ...and update the previous event */
	}
	
	x_ls = hhls(); /* compute LS solution */
	
	/* compute fitted values */
	for (i = 0; i < ESTW; i++) {
		b_hat[i] = 0.;
		for (j = 0; j < MNMC; j++)
			b_hat[i] += A_temp[i][j]*x_ls[j];
	}
	
	RSS = 0.; /* initialize */
	
	/* compute RSS */
	for (i = 0; i < ESTW; i++)
		RSS += (b_temp[i]-b_hat[i])*(b_temp[i]-b_hat[i]);
		
	return 1-RSS/TSS; /* return R^2 */

}

double calc_vrnc(long evnt) {
	
	short i, j; 
	double vrnc = 0., bhar_bar = 0., bhar_var = 0., temp = 0.;
	
	for (i = 0; i < ESTW; i++) {
		bhar_bar += b_temp[i]/ESTW;
		for (j = 0; j < MNMC; j++)
			bhar_bar -= A_temp[i][j]/(ESTW*MNMC);
	}
	
	for (i = 0; i < ESTW; i++) {
		temp = b_temp[i];
		for (j = 0; j < MNMC; j++)
			temp -= A_temp[i][j]/ESTW;
		bhar_var = (temp-bhar_bar)*(temp-bhar_bar)/ESTW;
	}
	
	return bhar_var;
}