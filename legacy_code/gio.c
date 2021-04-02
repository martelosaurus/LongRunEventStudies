/* titan input/output */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gio.h"
#include "glib.h"

extern M;
extern ret *R;
extern ymp *T;
extern double **P;
	
extern int ANMC;
extern int R_n;
extern int F_n;
extern int P_n;

int get_ret(char *file_name) {
	
	int year, mnth; 
	int rkey, firm;
	double rtrn, mrkt;
	
	FILE * fp;
	
	if ((fp = open_file(file_name,"r")) == NULL) 
		return 1;	
		
	fscanf(fp,"%*s");
		
	while(fscanf(fp,"%d,%d,%d,%d,%lf,%lf",&rkey,&firm,&year,&mnth,&rtrn,&mrkt) != EOF) {
	
		(R+R_n)->rkey = rkey;
		(R+R_n)->firm = firm;
		(R+R_n)->year = year;
		(R+R_n)->mnth = mnth;
		(R+R_n)->rtrn = rtrn;
		(R+R_n)->mrkt = mrkt;
	
		R_n++;
	}
	
	qsort(R,R_n,sizeof(ret),rkey_comp1);
	
	if (fclose(fp) != 0)
		return 1;
		
	return 0;
}

int get_ecf(char *file_name) {

	int excd, year, mnth;
	int i = 0, rkey, firm;
	double beta, size, btom;
	
	ecf *crnt;
	ymp key, *date;
	
	FILE * fp;

	if ((fp = open_file(file_name,"r")) == NULL) 
		return 1;
		
	fscanf(fp,"%*s");
		
	while(fscanf(fp,"%d,%d,%d,%d,%d,%lf,%lf",&rkey,&firm,&excd,&year,&mnth,&size,&btom) != EOF) {
	
		key.year = year; 
		key.mnth = mnth;
		
		if ((date = (ymp *) bsearch(&key,T,MNOM,sizeof(ymp),ymp_comp)) == NULL)
			continue;
		
		crnt = date->pntr+date->n;
		
		crnt->rkey = rkey;
		crnt->firm = firm;
		crnt->excd = excd;
		crnt->year = year;
		crnt->mnth = mnth;

		if (size < SIZE)
		    continue;
		
		if (calc_beta(&beta,crnt) == 1)
			continue;
	
		crnt->fctr  = (double *) malloc(MNMC*sizeof(double));
		crnt->rank  = (int *) malloc(MNMC*sizeof(int));

		if (crnt->fctr == NULL || crnt->rank == NULL)
			return 1;
		
		crnt->twin1 = (double *) malloc(MNMC*sizeof(double));
		crnt->twin2 = (double *) malloc(MNMC*sizeof(double));
		
		if (crnt->twin1 == NULL || crnt->twin2 == NULL)
			return 1;

		*(crnt->fctr+0) = beta;
		*(crnt->fctr+1) = size;
		*(crnt->fctr+2) = btom;
		
		if (++date->n == ECF)
			return 1;
	}
	
	if (fclose(fp) != 0)
		return 1;
	
	return 0;
}

int put_ecf(char *file_name) {
	
	int i, j;
	ecf *crnt;
	FILE *fp;
	
	if ((fp = open_file(file_name,"w")) == NULL)	
		return 1;
		
	fprintf(fp,"%s","rkey,fkey,firm,year,mnth,");
	fprintf(fp,"%s","beta,size,btom,");
	fprintf(fp,"%s","beta_rank,size_rank,btom_rank,");
	fprintf(fp,"%s","beta_twin,size_twin,btom_twin,");
	fprintf(fp,"%s",",,\n"); 
		
	for (i = 0; i < MNOM; i++)
		for (j = 0; j < (T+i)->n; j++) {
			crnt = (T+i)->pntr+j;
			fprintf(fp,"%d,%d,%d,%d,%d,",crnt->rkey,crnt->fkey,crnt->firm,crnt->year,crnt->mnth);
			fprintf(fp,"%f,%f,%f,",*(crnt->fctr+0),*(crnt->fctr+1),*(crnt->fctr+2));
			fprintf(fp,"%d,%d,%d,",*(crnt->rank+0),*(crnt->rank+1),*(crnt->rank+2));
			fprintf(fp,"%f,%f,%f,",*(crnt->twin1+0),*(crnt->twin1+1),*(crnt->twin1+2));
			fprintf(fp,"%f,%f,%f\n",*(crnt->twin2+0),*(crnt->twin2+1),*(crnt->twin2+2));
		}
	
	if (fclose(fp) != 0)
		return 1;
	
	return 0;
}

int put_pnl(char *file_name, int M_n) {

	int j;
	int i;
	FILE *fp;
	
	if ((fp = open_file(file_name, "w")) == NULL) {
		printf("can't open\n");
		return 1;
	}
	
	for (i = 0; i < P_n; i++) {
		for (j = 0; j < M_n+1; j++) 
			fprintf(fp,"%lf,",*(*(P+i)+j));
		fprintf(fp,"\n");
	}
		
	if (fclose(fp) != 0) {
		printf("can't close\n");
		return 1;
	}
	
	return 0;
}

FILE *open_file(char *file_name, char *mode) {
	
	char path[LEN];
	FILE *fp;
	strcpy(path,DIRECTORY);
	strcat(path,file_name);
	printf("path: %s\n", path);
	if ((fp = fopen(path,mode)) == NULL)
		return NULL;
	return fp;
}