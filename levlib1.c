#include <stdlib.h>
#include <math.h>

#include "levlib.h"

extern ymp *T;

int sort_fctr;
	
/* creates a hash table for the event/control firm panel */	
int new_ecf(void) {

	int i, j, k, F_n;
	ecf *F;
	htab ***H;
	
	for (i = 0; i < MNOM; i++) {
	
		F = (T+i)->pntr; 
		H = (T+i)->hash;
		F_n = (T+i)->n;
	
		/* base sort on rkey */
		qsort(F, F_n, sizeof(ecf), rkey_comp2);

		for (j = 0; j < F_n; j++)
			(F + j)->fkey = j;

		for (sort_fctr = 0; sort_fctr < MNMC; sort_fctr++) {
			qsort(F, F_n, sizeof(ecf), ecf_comp);
			for (k = 0; k < F_n; k++) {
				H[sort_fctr][0][k].a = H[sort_fctr][1][k].a = k;
				H[sort_fctr][0][k].b = H[sort_fctr][1][k].b = (F + k)->fkey;
			}
			qsort(H[sort_fctr][1], F_n, sizeof(htab), htab_compb);
		}
				
		/* base sort on fkey: do not change */
		qsort(F, F_n, sizeof(ecf), rkey_comp2);
		
		for (j = 0; j < MNMC; j++)
			if (calc_rank(T+i,j,QNTL) == 1)
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

int calc_rank(ymp *crnt, int fctr1, int q) {

	int X_n = 0;
	int i, j, F_n = crnt->n;
	double *qntl;
	double X[ECF];
	ecf *F = crnt->pntr;
	htab ***H = crnt->hash;

	if ((qntl = (double *)malloc(q*sizeof(double))) == NULL)
		return 1;

	if (!NYSE) {

		if (F_n >= q) {

			for (i = 0; i < q; i++) {
				j = (int)F_n*(i + 1.) / q - 1;
				qntl[i] = (F + H[fctr1][0][j].b)->fctr[fctr1];
			}

			for (i = 0; i < F_n; i++)
			for (j = q - 1; j >= 0; j--)
			if ((F + H[fctr1][0][i].b)->fctr[fctr1] <= qntl[j] && (F + H[fctr1][0][i].b)->rank[1] != -1)
				(F + H[fctr1][0][i].b)->rank[fctr1] = j;
		}
	}
	else {

		if (F_n >= q) {

			for (i = 0; i < F_n; i++)
			if ((F + H[fctr1][0][i].b)->excd == 1 && (F + H[fctr1][0][i].b)->rank[1] != -1) {
				X[X_n++] = (F + H[fctr1][0][i].b)->fctr[fctr1];
			}


			for (i = 0; i < q; i++) {
				j = (int)X_n*(i + 1.) / q - 1;
				qntl[i] = X[j];
			}

			for (i = 0; i < F_n; i++)
			for (j = q - 1; j >= 0; j--)
			if ((F + H[fctr1][0][i].b)->fctr[fctr1] <= qntl[j] && (F + H[fctr1][0][i].b)->rank[1] != -1)
				(F + H[fctr1][0][i].b)->rank[fctr1] = j;
		}
	}


	return 0;
}

int ecf_comp(const void *p1, const void *p2) {

	const double *a1 = ((const ecf *) p1)->fctr+sort_fctr;
	const double *a2 = ((const ecf *) p2)->fctr+sort_fctr;
	
	return (*a1>*a2)-(*a1<*a2);
}

int htab_compb(const void *p1, const void *p2) {

	const int *a1 = &(((const htab *)p1)->b);
	const int *a2 = &(((const htab *)p2)->b);

	return (*a1>*a2) - (*a1<*a2);
}