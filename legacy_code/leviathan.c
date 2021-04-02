#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "levio.h"
#include "levlib.h"
#include "perm.h"
#include "least_squares.h"

const mfs M[] = {
	/* {fctr1, fctr2, fctr2_upper, fctr2_lower, fctr3, fctr3_upper, fctr3_lower} */
	{ 0, -1, 1000., .5, -1, 0, 0 },
	{ 1, -1, 1000., .5, -1, 0, 0 },
	{ 2, -1, 1000., .5, -1, 0, 0 } 
};

const bss B = /* {fctr, factr_upper, fctr_lower} */ {-1, 0, 0};

#define ANMC 3
#define MAXR 1

ret *R;
ymp *T;	
double **P;
stack S;

int R_n = 0;
int F_n = 0;
int P_n = 0;

int driver(double *, ecf *, ymp *);
int new_stack(ecf *, ymp *);
int wrapper1(ecf *, grp *);
int wrapper2(ecf *, grp *); 

int main(void) {

	int i, j, k, flag;
	int start, finish; 
	ecf *evnt;
	
	printf("allocating memory for returns and event/control firm panels...\n");
	if (alloc1() == 1) {
		system("pause");
		return 1;
	}
	
	printf("reading returns panel from disk...\n");
	if (get_ret("grendel_ret.csv") == 1) {
		system("pause");
		return 1;
	}
	
	printf("reading event/control firm panel from disk...\n");
	if (get_ecf("grendel_ecf.csv") == 1) {
		system("pause");
		return 1;
	}
	
	printf("building event/control firm panel...\n");
	if (new_ecf() == 1) {
		system("pause");
		return 1;
	}
	
	printf("allocating memory for simulation panels...\n");
	if (alloc2(ANMC) == 1) {
		system("pause");
		return 1;
	}

	printf("simulating...\n");
	for (i = 0; i < MNOM; i++) {
		printf("year = %d: mnth = %d: P_n = %d:\n", (T + i)->year, (T + i)->mnth, P_n);
		for (j = 0; j < (T+i)->n; j++) {
		    evnt = (T+i)->pntr+j;
			
			if (evnt->fctr[0] == -99.)
				continue;
			
			if (B.fctr != -1 && (evnt->rank[B.fctr] > B.fctr_upper || evnt->rank[B.fctr] < B.fctr_lower))
				continue;
			if (driver(P[P_n],evnt,T+i) == 0)
				P_n++;
		}
	}

	printf("writing simulation panel to disk...\n");
	if (put_pnl("grendel_pnl.csv", ANMC) == 1) 
		return 1;

	if (ECF_OUT) { 
		printf("writing event/control firm panel to disk...\n");
		if (put_ecf("grendel_built_ecf.csv") == 1)
			return 1;
	}

	return 0;
}

int driver(double *bhar, ecf *crnt, ymp *date) {

	int i, j;

	int mnth = crnt->mnth%12+1, year = crnt->year+(mnth==1); /* next year-month */
	double evnt_rtrn[ANMC], ctrl_rtrn[ANMC]; /* event and control returns */
	ret *evnt[ANMC], *ctrl[ANMC]; /* pointers to returns panel for the event and control */
	grp *egrp, *cgrp; 
	
	for (i = 0; i < ANMC; i++) {
		evnt_rtrn[i] = 1.; 
		ctrl_rtrn[i] = 1.;
	}

	if (new_stack(crnt, date) == 1) /* load the stack */
		return 1;
	
	if ((egrp = stack_pop(&S)) == NULL || (cgrp = stack_pop(&S)) == NULL) 
		return 1;
		
	for (i = 0; i < FRWD; i++, mnth = mnth%12+1, year += (mnth==1)) {

		/* are all of the firms in the EVENT group active? */
		for (j = 0; j < ANMC; j++) {
			if ((evnt[j] = ret_srch(egrp->firm[j], i + 1, year, mnth)) == NULL) {
				j = -1;
				if ((egrp = stack_pop(&S)) == NULL)
					return 1;
			}
		}
			
		/* control */
		for (j = 0; j < ANMC; j++) {
			if ((ctrl[j] = ret_srch(cgrp->firm[j], i + 1, year, mnth)) == NULL) {
				j = -1;
				if ((cgrp = stack_pop(&S)) == NULL)
					return 1;
			}
		}
	
		/* returns */
		for (j = 0; j < ANMC; j++) {
			evnt_rtrn[j] *= (1.+evnt[j]->rtrn);
			ctrl_rtrn[j] *= (1.+ctrl[j]->rtrn);
		}
	}

	for (i = 0; i < ANMC; i++) {
		bhar[i] = evnt_rtrn[i]-ctrl_rtrn[i];
	}
	
	return 0;
}

int new_stack(ecf *crnt, ymp *date) {
	
	int i, j;
	static int **A = NULL;
	static ecf ***s = NULL;

	if (s == NULL) {

		if ((s = (ecf ***)malloc(ANMC*sizeof(ecf **))) == NULL)
			return 1;

		for (i = 0; i < ANMC; i++)
			if ((s[i] = (ecf **)malloc(MNCF*sizeof(ecf *))) == NULL)
				return 1;
	}
	
	S.N = perm_power(MNCF, ANMC)+1;

	if (A == NULL)
		A = perm(MNCF, ANMC);	
	
	for (i = 0; i < ANMC; i++)
		if (ecf_srch(s[i],crnt,date,M+i) == 1) /* load the firm stack */
			return 1;

	for (i = 0; i < S.N-1; i++) {
		for (j = 0; j < ANMC; j++)
			S.cgrp[i].firm[j] = s[j][A[i][j]];
		if (wrapper1(crnt, S.cgrp + i) == 1)
			return 1;
	}

	if (MAXR)
		stack_sort(&S);
	
	for (i = 0; i < ANMC; i++)
		S.cgrp[S.N - 1].firm[i] = crnt;
	
	S.n = S.N;

	return 0;
}

/* maximal R^2 */
int wrapper1(ecf *evnt, grp *cgrp) {

	int i, j, k, r;

	int mnth, year; /* month and year */
	const int crnt_mnth = 12 - (13 - evnt->mnth) % 12; /* last month */
	const int crnt_year = evnt->year - (crnt_mnth == 12); /* current (or last) year */
	double RSS = 0.; /* residual sum of squares */
	static double stock_bar, TSS;
	static double *stock = NULL, **portfolio = NULL; /* stock and portfolio returns */
	static double *stock_hat = NULL,  *stock_tmp = NULL, **portfolio_tmp = NULL; /* temporary stock and portfolio returns */
	ret *next; /* returns pointer */
	static ecf *last_stock, **last_portfolio = NULL; /* stock and portfolio trackers */

	/* allocate memory for the portfolio tracker */
	if (last_portfolio == NULL)
		if ((last_portfolio = (ecf **) malloc(ANMC*sizeof(ecf *))) == NULL)
			return 1;

	/* allocate memory for the stock returns */
	if (stock == NULL)
		if ((stock = (double *) malloc(BKWD*sizeof(double))) == NULL)
			return 1;

	/* allocate memory for the portfolio returns */
	if (portfolio == NULL) {

		if ((portfolio = (double **) malloc(BKWD*sizeof(double *))) == NULL)
			return 1;

		for (i = 0; i < BKWD; i++)
			if ((portfolio[i] = (double *) malloc((ANMC+1)*sizeof(double))) == NULL)
				return 1;

		for (i = 0; i < BKWD; i++)
			portfolio[i][ANMC] = 1.;
	}

	/* allocate memory for the predicted stock returns */
	if (stock_hat == NULL)
		if ((stock_hat = (double *) malloc(BKWD*sizeof(double))) == NULL)
			return 1;

	/* allocate memory for the temporary stock returns */
	if (stock_tmp == NULL)
		if ((stock_tmp = (double *) malloc(BKWD*sizeof(double))) == NULL)
			return 1;

	/* allocate memory for the temporary portfolio returns */
	if (portfolio_tmp == NULL) {

		if ((portfolio_tmp = (double **)malloc(BKWD*sizeof(double *))) == NULL)
			return 1;

		for (i = 0; i < BKWD; i++)
			if ((portfolio_tmp[i] = (double *)malloc((ANMC+1)*sizeof(double))) == NULL)
				return 1;

		for (i = 0; i < BKWD; i++)
			portfolio_tmp[i][ANMC] = 1.;
	}

	/* stock */
	if (evnt != last_stock) {

		/* time-step */
		i = 0; mnth = crnt_mnth; year = crnt_year;
		for (; i < BKWD; i++, mnth = 12 - (13 - mnth) % 12, year -= (mnth == 12)) {

			if ((next = ret_srch(evnt, -i - 1, year, mnth)) == NULL) {
				printf("returning on evnt\n");
				return 1;
			}

			stock[i] = next->rtrn;
		}

		stock_bar = 0.; 
		for (i = 0; i < BKWD; i++)
			stock_bar += stock[i] / BKWD;

		TSS = 0.;
		for (i = 0; i < BKWD; i++)
			TSS += (stock[i] - stock_bar) * (stock[i] - stock_bar);

		last_stock = evnt; /* update the stock tracker */
	}

	/* portfolio */
	for (i = 0; i < ANMC; i++) {

		if (cgrp->firm[i] != last_portfolio[i]) {

			/* time-step */
			j = 0; mnth = crnt_mnth; year = crnt_year;
			for (; j < BKWD; j++, mnth = 12 - (13 - mnth) % 12, year -= (mnth == 12)) {

				if ((next = ret_srch(cgrp->firm[i], -j - 1, year, mnth)) == NULL) {
					/*printf("returning on cgrp->firm[%d]\n", i);*/
					return 1;
				}

				portfolio[j][i] = next->rtrn; 
			}

			last_portfolio[i] = cgrp->firm[i]; /* update the portfolio tracker */
		}
	}
	
	/* bootstrap the stock and portfolio returns into the tmps */
	for (i = 0; i < BKWD; i++) {
		stock_tmp[i] = stock[i]; 
		for (j = 0; j < ANMC; j++)
			portfolio_tmp[i][j] = portfolio[i][j];
	}

	/* compute the betas and load them */
	if ((stock_tmp = hhls(portfolio_tmp, stock_tmp, BKWD, ANMC+1)) == NULL) 
		return 1;
	for (i = 0; i <= ANMC; i++) {
		cgrp->beta[i] = stock_tmp[i];
	}

	/* compute the predicted stock return */
	for (i = 0; i < BKWD; i++) {
		stock_hat[i] += cgrp->beta[ANMC];
		for (j = 0; j < ANMC; j++)
			stock_hat[i] += portfolio[i][j] * cgrp->beta[j];
	}

	/* compute the RSS*/
	for (i = 0; i < BKWD; i++)
		RSS += (stock_hat[i] - stock[i]) * (stock_hat[i] - stock[i]);

	cgrp->rsqr = 1 - RSS / TSS;

	return 0;
}

/* minimal skewness */
int wrapper2(ecf *evnt, grp *cgrp) {

	int i, j, k, r;

	int mnth, year; /* month and year */
	const int crnt_mnth = 12 - (13 - evnt->mnth) % 12; /* last month */
	const int crnt_year = evnt->year - (crnt_mnth == 12); /* current (or last) year */
	double mean, stdev, skew, hold1, hold2;
	double abnormal[BKWD];
	static double *stock = NULL, **portfolio = NULL; /* stock and portfolio returns */
	ret *next; /* returns pointer */
	static ecf *last_stock, **last_portfolio = NULL; /* stock and portfolio trackers */

	/* allocate memory for the portfolio tracker */
	if (last_portfolio == NULL)
		if ((last_portfolio = (ecf **)malloc(ANMC*sizeof(ecf *))) == NULL)
			return 1;

	/* allocate memory for the stock returns */
	if (stock == NULL)
		if ((stock = (double *)malloc(BKWD*sizeof(double))) == NULL)
			return 1;

	/* allocate memory for the portfolio returns */
	if (portfolio == NULL) {

		if ((portfolio = (double **)malloc(BKWD*sizeof(double *))) == NULL)
			return 1;

		for (i = 0; i < BKWD; i++)
			if ((portfolio[i] = (double *)malloc(ANMC*sizeof(double))) == NULL)
				return 1;
	}

	/* stock */
	if (evnt != last_stock) {

		/* time-step */
		i = 0; mnth = crnt_mnth; year = crnt_year;
		for (; i < BKWD; i++, mnth = 12 - (13 - mnth) % 12, year -= (mnth == 12)) {

			if ((next = ret_srch(evnt, -i - 1, year, mnth)) == NULL) 
				return 1;

			stock[i] = next->rtrn;
		}

		last_stock = evnt; /* update the stock tracker */
	}

	/* portfolio */
	for (i = 0; i < ANMC; i++) {

		if (cgrp->firm[i] != last_portfolio[i]) {

			/* time-step */
			j = 0; mnth = crnt_mnth; year = crnt_year;
			for (; j < BKWD; j++, mnth = 12 - (13 - mnth) % 12, year -= (mnth == 12)) {

				if ((next = ret_srch(cgrp->firm[i], -j - 1, year, mnth)) == NULL)
					return 1;

				portfolio[j][i] = next->rtrn;
			}

			last_portfolio[i] = cgrp->firm[i]; /* update the portfolio tracker */
		}
	}

	/* compute the abnormal returns */
	for (i = 0; i < BKWD; i++) {
		for (j = 0, hold1 = 1.; j < FRWD; j++) {
			r = (int) BKWD * (rand() - 1.) / RAND_MAX;
			hold1 *= (1. + stock[(r < 0) ? 0 : r]);
		}
		abnormal[i] = hold1;
		for (j = 0, hold1 = 1.; j < FRWD; j++) {
			r = (int) BKWD * (rand() - 1.) / RAND_MAX;
			for (k = 0, hold2 = 0.; k < ANMC; k++)
				hold2 += portfolio[r][j] / ANMC;
			hold1 *= (1. + hold2);
		}
	}

	/* compute the mean */
	for (i = 0, mean = 0.; i < BKWD; i++)
		mean += abnormal[i] / BKWD;

	/* compute the standard deviation */
	for (i = 0, stdev = 0.; i < BKWD; i++)
		stdev += (abnormal[i] - mean) * (abnormal[i] - mean) / (BKWD - 1);
	stdev = sqrt(stdev);

	/* compute the skewness */
	for (i = 0, skew = 0.; i < BKWD; i++)
		skew += (abnormal[i] - mean) * (abnormal[i] - mean) * (abnormal[i] - mean) / ((BKWD - 1) * stdev * stdev * stdev);

	cgrp->rsqr = 1. / skew;

	return 0;
}