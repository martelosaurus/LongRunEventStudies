#include <stdlib.h>
#include <math.h>

#include "glib.h"

extern ymp *T;

int sort_fctr;
	
/* creates a hash table for the event/control firm panel */	
int new_ecf(void) {

	int i, j, k, **H, F_n;
	ecf *F;
	
	for (i = 0; i < MNOM; i++) {
	
		F = (T+i)->pntr; 
		H = (T+i)->hash;
		F_n = (T+i)->n;
	
		/* base sort on rkey */
		qsort(F,F_n,sizeof(ecf),rkey_comp2);
	
		for (j = 0; j < F_n; j++)
			(F+j)->fkey = j;
		
		for (sort_fctr = 0; sort_fctr < MNMC; sort_fctr++) {
			qsort(F,F_n,sizeof(ecf),ecf_comp);
			for (k = 0; k < F_n; k++) 
				H[k][sort_fctr] = (F+k)->fkey;
		}
				
		/* base sort on fkey: do not change */
		qsort(F,F_n,sizeof(ecf),rkey_comp2);
		
		for (j = 0; j < MNMC; j++)
		    for (k = 1; k < F_n; k++)
				if (*((F+H[k-1][j])->fctr+j) > *((F+H[k][j])->fctr+j)) 
				    printf("warning!\n");
		
		for (j = 0; j < MNMC; j++)
			if (calc_rank(T+i,j) == 1)
				return 1;
	}
	
	return 0;
}

int calc_beta(double *beta, ecf *crnt) {

    static int warning = 1;
	int i, mnth = 12-(13-crnt->mnth)%12, year = crnt->year-(mnth==12); /* last year-month */
	double F[BKWD], M[BKWD], F_hat[BKWD];
	double F_bar = 0., M_bar = 0., cov_FF = 0., cov_FM = 0., cov_MM = 0.;
	ret *next;

	for (i = 0; i < BKWD; i++, mnth = 12-(13-mnth)%12, year -= (mnth==12)) {
		
		if ((next = ret_srch(crnt,-i-1,year,mnth)) == NULL)
			return 1;
		else {
			F[i] = next->rtrn;
				F_bar += F[i]/BKWD;
			M[i] = next->mrkt;
				M_bar += M[i]/BKWD;
		}
	}
	
	for (i = 0; i < BKWD; i++) {
		cov_FF += (F[i]-F_bar)*(F[i]-F_bar)/BKWD;
		cov_FM += (F[i]-F_bar)*(M[i]-M_bar)/BKWD;
		cov_MM += (M[i]-M_bar)*(M[i]-M_bar)/BKWD;
	}
	
	*beta = cov_FM / cov_MM;
	
	return 0;
}

int calc_rank(ymp *crnt, int fctr1) {

	int i, j, **H = crnt->hash, F_n = crnt->n;
	double qntl[QNTL]; 
	ecf *F = crnt->pntr;
			
	if (F_n >= QNTL) {

		for (i = 0; i < QNTL; i++) {
			j = (int) F_n*(i+1.)/QNTL-1;
			qntl[i] = *((F+H[j][fctr1])->fctr+fctr1);
		}
			
		for (i = 0; i < F_n; i++)
			for (j = QNTL-1; j >= 0; j--)
				if (*((F+H[i][fctr1])->fctr+fctr1) <= qntl[j])
					*((F+H[i][fctr1])->rank+fctr1) = j;
	}
		
	return 0;
}

int ecf_comp(const void *p1, const void *p2) {

	const double *a1 = ((const ecf *) p1)->fctr+sort_fctr;
	const double *a2 = ((const ecf *) p2)->fctr+sort_fctr;
	
	return (*a1>*a2)-(*a1<*a2);
}