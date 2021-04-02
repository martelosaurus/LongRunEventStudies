/* memory allocators */ 

#include <stdlib.h>
#include <math.h>

#include "glib.h"

extern ret *R;
extern ymp *T;
extern stack S;
extern double **P;

extern int R_n;

/* allocates memory for the returns and event/control firm panels */
int alloc1(void) {

	int i, j, **H, year = MINY, mnth = 1;
	
	if ((R = (ret *) malloc(RET*sizeof(ret))) == NULL)
		return 1;
		
	if ((T = (ymp *) malloc(MNOM*sizeof(ymp))) == NULL)
		return 1;
	
	for (i = 0; year <= MAXY; mnth = mnth%12+1, year += (mnth==1), i++) {
	
		(T+i)->year = year;
		(T+i)->mnth = mnth;
		
		if (((T+i)->pntr = (ecf *) malloc(ECF*sizeof(ecf))) == NULL)
			return 1;
			
		if (((T+i)->hash = (int **) malloc(ECF*sizeof(int *))) == NULL)
			return 1;
		
		for (j = 0; j < ECF; j++) 
			if ((*((T+i)->hash+j) = (int *) malloc(MNMC*sizeof(int))) == NULL)
				return 1;
				
		(T+i)->n = 0;
	}
	
	return 0;
}

/* allocate stack AND it's groups */

/* allocates memory for the BHAR panel */
int alloc2(int M_n) {

	int i;
	
	if ((P = (double **) malloc(ECF*MNOM*sizeof(double *))) == NULL)
		return 1;
		
	for (i = 0; i < ECF*MNOM; i++)
		if ((*(P+i) = (double *) malloc(M_n*sizeof(double))) == NULL)
			return 1;
		
	S.n = 0; 
	S.N = perm_power(MNCF, M_n)+1;

	if ((S.cgrp = (grp *)malloc(S.N*sizeof(grp))) == NULL)
		return 1;

	for (i = 0; i < S.N; i++)
		if ((S.cgrp[i].firm = (ecf **) malloc(M_n*sizeof(ecf *))) == NULL)
			return 1;
		
	return 0;
}