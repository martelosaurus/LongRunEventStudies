/* ECFP */

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
	 
extern ecf *T; /* ECFP... */
	extern long T_n; /* #of observations */

	extern ecf *T_size; /* ...sorted by size */
	extern ecf *T_btom; /* ...sorted by btom */
	extern ecf *T_beta; /* ...sorted by beta */
	extern ecf *T_sicc; /* ...sorted by sicc */
	extern ecf *T_mome; /* ...sorted by mome */
	
extern ecd *D; /* event/control firm deciles */
	extern short D_n; /* counter */

extern long **C; /* controls: C: MNCF x MNMC */
	extern short C_n; /* control counter */

/* creates sorted ECFPs and computes decile ranks */
void ecf_build(void) {
	
	long i; /* counter */
	
	/* create key for ECFP sorted by ret_key: T_n < R_n so glo_key != ret_key */
	for (i = 0; i < T_n; i++) {
		(T+i)->glo_key = i;
		(T+i)->used = -1;
		(T+i)->size_twin = -1;
		(T+i)->btom_twin = -1;
		(T+i)->beta_twin = -1;
	}

	/* allocate memory for sorted ECFPs */
	T_size = (ecf *) malloc(T_n*sizeof(ecf));
	T_btom = (ecf *) malloc(T_n*sizeof(ecf));
	T_beta = (ecf *) malloc(T_n*sizeof(ecf));
	T_sicc = (ecf *) malloc(T_n*sizeof(ecf));
	T_mome = (ecf *) malloc(T_n*sizeof(ecf));
	
	if (T_size == NULL || T_btom == NULL || T_beta == NULL || T_sicc == NULL || T_mome == NULL)
		general_error("error: mem_all: ecf_sort");

	/* create sorted copies of ECFPs */
	ecf_sort(T_size, size_comp);
	ecf_sort(T_btom, btom_comp); 
	ecf_sort(T_beta, beta_comp); 
	ecf_sort(T_sicc, sicc_comp); 
	ecf_sort(T_mome, mome_comp); 

	/* allocate memory for decile object */
	if ((D = (ecd *) malloc(10*MNOM*sizeof(ecd))) == NULL)
		general_error("error: mem_all: ecf_build"); 
		
	/* compute deciles */
	ecf_decls(T_size,1,size_decl_pick);
	ecf_decls(T_btom,1,btom_decl_pick);
	ecf_decls(T_mome,1,mome_decl_pick);
	
	/* assign decile ranks */
	ecf_ranks(T);		
	ecf_ranks(T_size);
	ecf_ranks(T_btom);
	ecf_ranks(T_beta);
	ecf_ranks(T_sicc);
	ecf_ranks(T_mome);

}

/* copies observations from 'from' to 'to' */
long ecf_copy(ecf *to, ecf *from, short nyse, short year, short mnth) {
	
	long i, to_n = 0; /* counter and 'to' counter */
	
	for (i = 0; i < T_n; i++) {
	
		/* if filtering on exchange and exchange not NYSE, continue */
		if (nyse && from[i].excd != 1)
			continue;
		
		/* if filtering on year-month and year-month incorrect, continue */
		if (year>0 && mnth>0 && year != from[i].year && mnth != from[i].mnth)
			continue;
			
		/* copy... */ 	
		to[to_n].ret_key = from[i].ret_key;
		to[to_n].glo_key = from[i].glo_key;
		to[to_n].used = from[i].used;
		to[to_n].firm = from[i].firm;
		to[to_n].excd = from[i].excd;
		to[to_n].year = from[i].year;
		to[to_n].mnth = from[i].mnth;
		to[to_n].size = from[i].size;
		to[to_n].btom = from[i].btom;
		to[to_n].beta = from[i].beta;
		to[to_n].sicc = from[i].sicc;
		to[to_n].mome = from[i].mome;
	
		to_n++; /* ...and update the counter */
	}
	
	return to_n; /* return the # copied observations */
}

/* creates sorted copies of ECFP */
void ecf_sort(ecf *T_mach, int (*comp)(const void *, const void *)) {

	/* np: new ECFP */

	long i; /* counter */
	
	/* copy contents of base ECFP to new ECFP */
	ecf_copy(T_mach,T,0,0,0);
	
	/* sort new ECFP by its factor */
	qsort(T_mach,T_n,sizeof(ecf),comp);
		
	/* create key for new ECFP */
	for (i = 0; i < T_n; i++)
		(T_mach+i)->loc_key = i;
}

/* computes deciles for size, btom and mome */
void ecf_decls(ecf *T_mach, short nyse, void (*decl_pick) (ecd *, const ecf *)) {

	short i, rank, year = MINY-1, mnth = 0; /* counter, rank, year and month */
	long T_sub_n; /* #of sub-sample ECFP observations */

	ecf *T_sub; /* sub-sample ECFP */
	
	/* allocate memory for the sub-sample ECFP */
	if ((T_sub = (ecf *) malloc(MECF*sizeof(ecf))) == NULL)
		general_error("error: mem_all: ecf_decls");

	/* loop over the months in the sample window */
	for (i = 0; i < MNOM; i++) {
		
		/* update year and mnth */
		if ((mnth = mnth%12+1) == 1)
			year++;	
		
		/* copy the sub-sample ECFP observations and count them */
		T_sub_n = ecf_copy(T_sub,T_mach,nyse,year,mnth);
		
		/* loop over the decile ranks */
		for (rank = 0; rank < 10; rank++) {
			(D+rank+10*i)->year = year; /* update year */
			(D+rank+10*i)->mnth = mnth; /* update month */
			(D+rank+10*i)->rank = rank; /* update rank */ 
			decl_pick(D+rank+10*i,T_sub+(long) (.1*(rank+1.)*T_sub_n+.5)-1); /* pick the percentile ECFP */
		}
	}
	
	free(T_sub); /* free the sub-sample ECPF */
}

/* computes decile ranks for each observation in ECFP */
void ecf_ranks(ecf *T_mach) {
	
	short i, rank, year = MINY-1, mnth = 0; /* counter, rank, year and month */
	long j; /* counter */
	
	/* loop over the months in the sample window */
	for (i = 0; i < MNOM; i++) {
	
		/* update year and mnth */
		if ((mnth = mnth%12+1) == 1)
			year++;	
	
		/* loop over the ECFP observations */
		for (j = 0; j < T_n; j++) {
		
			/* if incorrect year or month, continue */
			if ((T_mach+j)->year != year || (T_mach+j)->mnth != mnth)
				continue;
		
			/* assign each matching char. the correct rank */
			for (rank = 9; rank > -1; rank--) {				
				if ((T_mach+j)->size <= D[rank+10*i].size) 
					(T_mach+j)->size_rank = rank;
				if ((T_mach+j)->btom <= D[rank+10*i].btom)
					(T_mach+j)->btom_rank = rank;
				if ((T_mach+j)->mome <= D[rank+10*i].mome)
					(T_mach+j)->mome_rank = rank;
			}	
		}
	}
}

/* searches decile object for desired decile */
double decl_srch(char mach, short rank, short year, short mnth) {
	int i; /* counter */
	for (i = 0; i < 10*MNOM; i++)
		if ((D+i)->year == year && (D+i)->mnth == mnth && (D+i)->rank == rank) {
			switch(mach) {
				case 's': return (D+i)->size; /* return the size decile */
				case 'b': return (D+i)->btom; /* return the btom decile */ 
				case 'm': return (D+i)->mome; /* return the mome decile */
			}
	}
	return -1.;
}
		
/* searches the specified ECFP for a set of controls for the specified event */ 
long ecf_srch(ecf *T_mach, long evnt, char decl_mach, void (* pick)(ecf *, const ecf *), int (* comp)(const void *, const void *)) {

	/*	T_mach: a sorted ECFP 
		evnt: GLO_KEY for current event
		decl_mach: decile matching characteristic
		PICK: pointer to pick function
		COMP: pointer to comparison function
	*/
	
	long lsrv; 
	short i; /* counter */
	
	ecf key, *ptr; /* search key and pointer to result */

	/* copy factor value of current event into the search key */
	pick(&key,T+evnt); 

	/* search for the factor value in the factor sorted ECFP... */
	if ((ptr = (ecf *) bsearch(&key,T_mach,T_n,sizeof(ecf),comp)) == NULL)
		return -1; /* ...and return an error code if it can't be found (shouldn't happen) */

	/* use the result of the above search to search for suitable controls... */
	if ((lsrv = ecf_loc_srch(T_mach,(ptr->loc_key),evnt,decl_mach)) == -1)
		return -1; /* ...and return an error code if an insufficient number of controls were found */
	
	/* otherwise... */
	C_n++; /* ...update the matching char. counter... */
	return lsrv; /* ...and return 0 */
}

/* subroutine of ecf_srch: performs local-linear search */
long ecf_loc_srch(ecf *T_mach, long evnt_loc, long evnt_glo, char decl_mach) {
 
	/*	evnt_loc: LOC_KEY of current event 
		evnt_glo: GLO_KEY of current event
	*/
 
	long i, k, first = -1; /* search counter */
	short C_m = 0, flag, s; /* controls counter, counter and sign */
	short same_year, same_mnth, diff_firm, same_rank, wsnt_used = 1; /* flags */

	/* this loop progresses like -1,1,-2,2,-3,3,... and stops if the local key es+s*k is out of bounds */
	for (k = 0, s = -1; C_m < MNCF && k < T_n; k += (s+1)/2, s *= -1) {
	
		if (evnt_loc+s*k < T_n && evnt_loc+s*k >= 0) {
		
			//if ((wsnt_used = ((T_mach+evnt_loc+s*k)->used == -1)) == 0) continue; /* wasn't used? */  
			if ((same_year = ((T_mach+evnt_loc+s*k)->year == (T+evnt_glo)->year)) == 0) continue; /* same year? */
			if ((same_mnth = ((T_mach+evnt_loc+s*k)->mnth == (T+evnt_glo)->mnth)) == 0) continue; /* same mnth? */
			if ((diff_firm = ((T_mach+evnt_loc+s*k)->firm != (T+evnt_glo)->firm)) == 0) continue; /* different firm? */

			/* set the same_rank flag value equal to 1 if the candidate is in the same decile */
			switch(decl_mach) {
				case 's': 
					same_rank = ((T_mach+evnt_loc+s*k)->size_rank == (T+evnt_glo)->size_rank);
					break;
				case 'b': 
					same_rank = ((T_mach+evnt_loc+s*k)->btom_rank == (T+evnt_glo)->btom_rank);
					break;
				case 'm':
					same_rank = ((T_mach+evnt_loc+s*k)->mome_rank == (T+evnt_glo)->mome_rank);
					break;
				case 'l': /* for Lyon, Barber and Tsai */
					same_rank = (.7*(T+evnt_glo)->size < (T_mach+evnt_loc+s*k)->size && (T_mach+evnt_loc+s*k)->size < 1.3*(T+evnt_glo)->size);
					break;
				default: 
					same_rank = 1;
					break;
			}
			
			/* if a suitable control has been found, store its return key */
			if (wsnt_used && same_year && same_mnth && diff_firm && same_rank) {
			
				*(*(C+C_m)+C_n) = (T+(T_mach+evnt_loc+s*k)->glo_key)->ret_key; /* store return key */
				//(T_mach+evnt_loc+s*k)->used = 1;
				
				if (first == -1)
					first = evnt_loc+s*k;
				
				C_m++; /* otherwise, update the control count */
			}
		}
	}
	
	/* return -1 if an insufficient number of controls were found...
	   and the number of  otherwise */
	return (C_m<MNCF)?-1:first;
}