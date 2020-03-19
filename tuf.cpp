/* titan utilities and other functions */

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

/* computes the sample covariance of x and y */
/* not intended for use outside of beta calculation */
double cov(double *x, double *y) {

	long i; /* counter */
	double s_xy, x_bar, y_bar; /* sample covariance and means */ 
	s_xy = x_bar = y_bar = 0.; /* initialize */
		
	/* compute the sample means */
	for (i = 0; i < ESTW; i++) {
		x_bar += *(x+i)/ESTW;
		y_bar += *(y+i)/ESTW;
	}
	
	/* compute the sample covariance */
	for (i = 0; i < ESTW; i++)
		s_xy += (*(x+i)-x_bar)*(*(y+i)-y_bar)/ESTW;
		
	/* return the sample covariance */
	return s_xy;
}

/* reports general error */
void general_error(char *error) {
	
	printf("%s\n", error); 
	system("pause"); 
	exit(EXIT_FAILURE); 
}

/* reports error in opening a file */
void open_error(char *file_name) {
	
	char error[LLEN];
	strcpy(error, "error: unable to open ");
	strcat(error, file_name);
	printf("%s\n", error);
	system("pause");
	exit(EXIT_FAILURE);
}

/* reports error in closing a file */
void close_error(char *file_name) {
	
	char error[LLEN];
	strcpy(error, "error: unable to close ");
	strcat(error, file_name);
	printf("%s\n", error);
	system("pause");
	exit(EXIT_FAILURE);
}

/* allocates memory for an m x 1 array of longs */
short *s1alloc(long m) {
	short *p;
	if ((p = (short *) malloc(m*sizeof(short))) == NULL)
		general_error("error: mem_all: l1alloc");
	return p;
}

/* allocates memory for an m x 1 array of longs */
long *l1alloc(long m) {
	long *p;
	if ((p = (long *) malloc(m*sizeof(long))) == NULL)
		general_error("error: mem_all: l1alloc");
	return p;
}

/* allocates memory for an m x 1 array of floats */
double *f1alloc(long m) {
	double *p; 
	if ((p = (double *) malloc(m*sizeof(double))) == NULL)
		general_error("error: mem_all: f1alloc");
	return p;
}

/* allocates memory for an m x n array of longs */
long **l2alloc(long m, long n) {
	long i, **p;
	if ((p = (long **) malloc(m*sizeof(long *))) == NULL)
		general_error("error: mem_all: lalloc: **p");
	for (i = 0; i < m; i++)
		*(p+i) = l1alloc(n);
	return p;
}

/* allocates memory for an m x n array of floats */
double **f2alloc(long m, long n) {
	long i;
	double **p;
	if ((p = (double **) malloc(m*sizeof(double *))) == NULL)
		general_error("error: mem_all: lalloc: **p");
	for (i = 0; i < m; i++)
		*(p+i) = f1alloc(n);
	return p;
}

short rank_pick(ecf *x, char mach) {
	switch(mach) {
		case 's': return x->size;
		case 'b': return x->btom;
		case 'm': return x->mome;
		default: return -1;
	}
}

/* assigns size value from 'from' to 'to' */
void size_decl_pick(ecd *to, const ecf *from) {
	to->size = from->size;
}

/* assigns size value from 'from' to 'to' */
void btom_decl_pick(ecd *to, const ecf *from) {
	to->btom = from->btom;
}

/* assigns size value from 'from' to 'to' */
void mome_decl_pick(ecd *to, const ecf *from) {
	to->mome = from->mome;
}

/* assigns size value from 'from' to 'to' */
void size_pick(ecf *to, const ecf *from) {
	to->size = from->size;
}

/* assigns btom value from 'from' to 'to' */
void btom_pick(ecf *to, const ecf *from) {
	to->btom = from->btom;
}

/* assigns beta value from 'from' to 'to' */
void beta_pick(ecf *to, const ecf *from) {
	to->beta = from->beta;
}

/* assigns sicc value from 'from' to 'to' */
void sicc_pick(ecf *to, const ecf *from) {
	to->sicc = from->sicc;
}

/* assigns mome value from 'from' to 'to' */
void mome_pick(ecf *to, const ecf *from) {
	to->mome = from->mome;
}

/* firm comparison function for ECFP */
int ecf_ret_key_comp(const void *p1, const void *p2) {

	const long *a1 = &(((const ecf *)p1)->ret_key);
	const long *a2 = &(((const ecf *)p2)->ret_key);
	
	return (*a1>*a2)-(*a1<*a2);
}

/* firm comparison function for returns panel */
int ret_ret_key_comp(const void *p1, const void *p2) {

	const long *a1 = &(((const ret *)p1)->ret_key);
	const long *a2 = &(((const ret *)p2)->ret_key);
	
	return (*a1>*a2)-(*a1<*a2);
}

/* size comparison function for ECFP */
int size_comp(const void *p1, const void *p2) {

	const double *a1 = &(((const ecf *)p1)->size);
	const double *a2 = &(((const ecf *)p2)->size);
	
	return (*a1>*a2)-(*a1<*a2);
}

/* btom comparison function for ECFP */
int btom_comp(const void *p1, const void *p2) {

	const double *a1 = &(((const ecf *)p1)->btom);
	const double *a2 = &(((const ecf *)p2)->btom);
	
	return (*a1>*a2)-(*a1<*a2);
}

/* beta comparison function for ECFP */
int beta_comp(const void *p1, const void *p2) {

	const double *a1 = &(((const ecf *)p1)->beta);
	const double *a2 = &(((const ecf *)p2)->beta);
	
	return (*a1>*a2)-(*a1<*a2);
}

/* sicc comparison function for ECFP */
int sicc_comp(const void *p1, const void *p2)  {

	const short *a1 = &(((const ecf *)p1)->sicc);
	const short *a2 = &(((const ecf *)p2)->sicc);
	
	return (*a1>*a2)-(*a1<*a2);
}	

/* mome comparison function for ECFP */
int mome_comp(const void *p1, const void *p2)  {

	const double *a1 = &(((const ecf *)p1)->mome);
	const double *a2 = &(((const ecf *)p2)->mome);
	
	return (*a1>*a2)-(*a1<*a2);
}	