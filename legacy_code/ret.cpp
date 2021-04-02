/* returns panel */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ecf.h" 
#include "ret.h" 
#include "tio.h"
#include "tls.h" 
#include "tps.h" 
#include "ttc.h"
#include "tuf.h"

extern ret * R; /* returns panel */
	extern long R_n; /* # observations */

	extern ret_mkt *R_mkt; /* market returns panel */
		extern short R_mkt_n; /* # observations */
		
	extern ret_ref *R_size; /* RPRP */
	extern ret_ref *R_btom;
	extern ret_ref *R_mome; 
		extern short R_ref_n; /* # observations */
		
extern ecf *T; /* ECFP... */
	extern long T_n; /* #of observations */

/* computes an event/control's beta */
short calc_beta(double *beta, long evnt) {

	/* evnt: RET_KEY for current event */ 

	short i, year = (R+evnt)->year, mnth = (R+evnt)->mnth; /* counter, year and month */
	long firm = (R+evnt)->firm;
	double *frets, *mrets; /* firm returns and market returns */
	
	/* allocate memory for firm returns and market returns */
	frets = f1alloc(ESTW); 
	mrets = f1alloc(ESTW);
	
	/* loop over the months in the holding period */
	for (i = 0; i < ESTW; i++) {
	
		/* update year and mnth */
		if ((mnth = 12-(13-mnth)%12) == 12)
			year--;
			
		/* search for the event/control in the returns panel... */
		if (ret_srch(frets+i,evnt-i-1,year,mnth,firm) == 1) 
			return 1; /* ...and return an error code if it can't be found */
	
		/* search for the event/control in the returns panel... */
		if (ret_mkt_srch(mrets+i,year,mnth) == 1) 
			return 1; /* ...and return an error code if it can't be found */
	}

	/* compute beta */
	*beta = cov(frets,mrets)/cov(mrets,mrets); 

	free(frets); free(mrets); /* free the firm returns and market returns */
	
	/* all ok */
	return 0;
}

/* computes an event/control's momentum */
short calc_mome(double *mome, long evnt) {
	
	/* evnt: RET_KEY for current event */ 
	
	short i, year = (R+evnt)->year, mnth = (R+evnt)->mnth; /* counter, year and month */
	long firm = (R+evnt)->firm; /* firm */
	double fret; /* firm return */
	
	*mome = 1.; /* initialize */
	
	/* loop over the months in the holding period */
	for (i = 0; i < MOMW; i++) {
	
		/* update year and mnth */
		if ((mnth = 12-(13-mnth)%12) == 12)
			year--;
		
		/* search for the event/control in the returns panel... */
		if (ret_srch(&fret,evnt-i-1,year,mnth,firm) == 1)
			return 1; /* ...and return an error code if it can't be found */
		
		/* update the momentum */
		*mome *= (1.+fret);
	}
	
	/* all ok */
	return 0;
}

/*
void ref_build(void) {
	
	R_size = (ret_ref *) malloc(10*MNOM*sizeof(ret_ref));
	R_btom = (ref_ref *) malloc(10*MNOM*sizeof(ret_ref));
	R_mome = (ret_ref *) malloc(10*MNOM*sizeof(ret_ref));
	
	if (R_size == NULL || R_btom == NULL || R_mome == NULL)
		general_error("error: mem_all: ref_build");

	ref_calc(R_size);
	ref_calc(R_btom);
	ref_calc(R_mome);
}

void ref_calc(ret_ref *R_ref) {

	short i, j, year, mnth; 
	long k, pret_n; 
	double pret;
	
	for (i = 0; i < 10; i++) { 
			
		year = MINY; mnth = 0; 
			
		for (j = 0; j < MNOM; j++) { 
				
			pret = 0.; 
				
			if ((mnth = mnth%12+1) == 1)
				year = year+1;
					
			pret_n = 0;

			for (k = 0; k < T_n; k++) 
				if (rank_pick(T+k,'s') == i) {
					pret += (R+(T+k)->ret_key)->fret;
					pret_n++; 
				}
					
			(R_ref+R_ref_n)->rank = i;
			(R_ref+R_ref_n)->year = year;
			(R_ref+R_ref_n)->mnth = mnth;
			(R_ref+R_ref_n)->pret = pret/pret_n; 
					
			R_ref_n++;			
		}
	}
}
/*

/* searches the returns panel for the specified event/control and updates firm returns */
short ret_srch(double *fret, long ret_key, short year, short mnth, long firm) {

	/* ret_key: starting ret_key (not necessarily event's ret_key) */ 
 
	long k, s; /* search counter, sign, firm and return key */
	short crct_year, crct_mnth, crct_firm; /* flags */
	
	/* this loop progresses like -1,1,-2,2,-3,3,... and stops if the local key es+s*k is out of bounds */
	for (k = 0, s = -1; k < MNOM; k += (s+1)/2, s *= -1) {
	
		if (ret_key+s*k < R_n && ret_key+s*k >= 0) {
	
			crct_year = ((R+ret_key+s*k)->year == year); /* correct year? */
			crct_mnth = ((R+ret_key+s*k)->mnth == mnth); /* correct mnth? */
			crct_firm = ((R+ret_key+s*k)->firm == firm); /* correct firm? */
		 
			/* if the desired return has been found, return a pointer to it */
			if (crct_year && crct_mnth && crct_firm) {
				*fret = (R+ret_key+s*k)->fret;
				return 0;
			}
		}
	}
	
	/* otherwise, return an error code */
	return 1;
}

short ret_mkt_srch(double *mret, short year, short mnth) {
	short i;
	for (i = 0; i < MNOM; i++)
		if ((R_mkt+i)->year == year && (R_mkt+i)->mnth == mnth) {
			*mret = (R_mkt+i)->mret;
			return 0;
		}
	return 1;
}

/*
short ret_ref_srch(char type, double *pret, short rank, short year, short mnth) {	
	long i;
	for (i = 0; i < MNRP*MNOM; i++)
		if ((R_ref+i)->size_rank == size_rank && (R_ref+i)->year == year && (R_ref+i)->mnth == mnth) {
			
		}
		
	return 1;
}
*/