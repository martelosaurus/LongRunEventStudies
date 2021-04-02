/* titan input/output */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ecf.h" 
#include "ret.h" 
#include "tio.h"
#include "tls.h" 
#include "tps.h" 
#include "ttc.h"
#include "tuf.h"

extern ret *R; /* returns panel */
	extern long R_n; /* # observations */
		 
extern ret_mkt *R_mkt; /* market returns panel */
	extern short R_mkt_n; /* # observations */
		 
extern ecf *T; /* ECFP... */
	extern long T_n; /* # observations */

extern double **P; /* bhar panel */
	extern long P_n; /* # observations */		
	
extern long *Q; /* event dates */
	
/* reads the ECFP from disk */
void ecf_read(char file_name[]) {

	long i = 0; /* counter */

	char path[LLEN],line[LLEN]; /* file path and scanning string */
	
	/* input variables */
	short mnth,year,sicc,excd;
	long ret_key,firm;
	double size,btom,beta,mome;
	
	FILE * fp; /* file pointer */

	/* allocate memory for the ECFP */
	if ((T = (ecf *) malloc(MECF*sizeof(ecf))) == NULL) 
		general_error("error: mem_all: ecf_read");  
	
	/* create full file name */
	strcpy(path,DIRECTORY); 
	strcat(path,file_name);
	
	/* open file */
	if ((fp = fopen(path,"r")) == NULL) 
		open_error(file_name);
	
	/* scan contents of file into scanning string */
	while(fgets(line,LLEN,fp) != NULL) {
	
		/* parse the scanning string */
		if (sscanf(line,"%ld,%ld,%hd,%hd,%hd,%lf,%lf,%hd",&ret_key,&firm,&excd,&year,&mnth,&size,&btom,&sicc) == 8) {	
			
			/* compute beta */
			if (calc_beta(&beta,ret_key) == 1)
				continue;
			
			/* compute mome */
			if (calc_mome(&mome,ret_key) == 1)
				continue;
		
			/* load ECFP */
			(T+T_n)->ret_key = ret_key;
			(T+T_n)->firm = firm;
			(T+T_n)->excd = excd;
			(T+T_n)->mnth = mnth;
			(T+T_n)->year = year;
    		(T+T_n)->size = size;
			(T+T_n)->btom = btom;
			(T+T_n)->beta = beta;
			(T+T_n)->sicc = sicc;
			(T+T_n)->mome = mome;
			
			T_n++; 
        }
	}
	
	/* sort ECFP by firm */
	qsort(T,T_n,sizeof(ecf),ecf_ret_key_comp);
	
	/* close file */
	if (fclose(fp) != 0)
		close_error(file_name);
}

/* writes the ECFP to disk */
void ecf_write(char file_name[]) {

	long i = 0; /* counter */

	char path[LLEN],line[LLEN]; /* file path and scanning string */
	char head[] = "glo_key,firm,year,mnth,size,btom,beta,sicc,mome,size_rank,btom_rank,mome_rank,size_twin,btom_twin,beta_twin\n"; /* header */
	char frmt[] = "%ld,%ld,%hd,%hd,%10.16f,%10.15f,%10.15f,%hd,%10.15f,%hd,%hd,%hd,%ld,%ld,%ld\n"; /* format */
	
	FILE * fp; /* file pointer */
	
	/* create full file name */
	strcpy(path,DIRECTORY); 
	strcat(path,file_name);
	
	/* open file */
	if ((fp = fopen(path,"w")) == NULL) 
		open_error(file_name);
	
	/* write each line of the ECFP to the file */
	fprintf(fp,head);
	for (i = 0; i < T_n; i++) 
		fprintf(fp,frmt,T[i].glo_key,T[i].firm,T[i].year,T[i].mnth,T[i].size,T[i].btom,T[i].beta,T[i].sicc,T[i].mome,T[i].size_rank,T[i].btom_rank,T[i].mome_rank,T[i].size_twin,T[i].btom_twin,T[i].beta);
	
	/* close file */
	if (fclose(fp) != 0)
		close_error(file_name);
}

/* reads the returns panel from disk */
void ret_read(char file_name[]) {

	long i; /* counter */
	char path[LLEN],line[LLEN]; /* file path and scanning string */
	
	/* input variables */
	short mnth,year; 
	long firm,ret_key;
	double fret,mret;
	
	FILE * fp; /* file pointer */

	/* allocate memory for the returns panel */
	if ((R = (ret *) malloc(MRET*sizeof(ret))) == NULL)
		general_error("error: mem_all: ret_read");  
	
	/* create full file name */
	strcpy(path,DIRECTORY); 
	strcat(path,file_name);
	
	/* open file */
	if ((fp = fopen(path,"r")) == NULL) 
		open_error(file_name);
	
	/* scan contents of file into scanning string */
	while(fgets(line,LLEN,fp) != NULL) {
	
		/* parse the scanning string */
		if (sscanf(line,"%ld,%ld,%hd,%hd,%lf",&ret_key,&firm,&year,&mnth,&fret) == 5) {	
		
			/* load returns panel */
			(R+R_n)->ret_key = ret_key;
			(R+R_n)->firm = firm;
			(R+R_n)->mnth = mnth;
			(R+R_n)->year = year;
			(R+R_n)->fret = fret;
			
			R_n++;
        }
	}
	
	/* sort returns panel by firm */
	qsort(R,R_n,sizeof(ret),ret_ret_key_comp);

	/* check that the ret_keys have loaded properly */
	for (i = 1; i < R_n; i++)
		if (i != R[i].ret_key) 
			general_error("error: : ret_read");

	/* close file */
	if (fclose(fp) != 0)
		close_error(file_name);
}

/* reads the market returns panel from disk */
void ret_mkt_read(char file_name[]) {

	long i; /* counter */
	char path[LLEN],line[LLEN]; /* file path and scanning string */
	
	/* input variables */
	short year, mnth;
	double mret;
	
	FILE * fp; /* file pointer */

	/* allocate memory for the returns panel */
	if ((R_mkt = (ret_mkt *) malloc(MNOM*sizeof(ret))) == NULL)
		general_error("error: mem_all: ret_read");  
	
	/* create full file name */
	strcpy(path,DIRECTORY); 
	strcat(path,file_name);
	
	/* open file */
	if ((fp = fopen(path,"r")) == NULL) 
		open_error(file_name);
	
	/* scan contents of file into scanning string */
	while(fgets(line,LLEN,fp) != NULL) {
	
		/* parse the scanning string */
		if (sscanf(line,"%hd,%hd,%lf",&year,&mnth,&mret) == 3) {	
		
			/* load returns panel */
			(R_mkt+R_mkt_n)->year = year;
			(R_mkt+R_mkt_n)->mnth = mnth;
			(R_mkt+R_mkt_n)->mret = mret; 
			
			R_mkt_n++;
        }
	}

	/* close file */
	if (fclose(fp) != 0)
		close_error(file_name);
}

/* writes the BHAR panel to disk */
void pnl_write(char file_name[]) {

	short j; /* counter */
	long i; /* counter */ 

	char path[LLEN]; /* file path and scanning string */
	
	/* create full file name */
	strcpy(path,DIRECTORY); 
	strcat(path,file_name);
	
	FILE * fp; /* file pointer */
	
	/* open file */
	if ((fp = fopen(path,"w")) == NULL) 
		open_error(file_name);
	
	/* loop over the BHARs */
	for (i = 0; i < P_n; i++) {
		
		j = 0; /* reset the matching char. counter */
		
		if (i < 10) {
			printf(BPOF,**(P+i));
			printf("\n");
		}
		
		if (MACH || MRSQ) {
			if (SIZE) fprintf(fp,BPOF,*(*(P+i)+j++));
				if (SIZE_BTOM) fprintf(fp,BPOF,*(*(P+i)+j++));
				if (SIZE_MOME) fprintf(fp,BPOF,*(*(P+i)+j++));
			if (BTOM) fprintf(fp,BPOF,*(*(P+i)+j++));
				if (BTOM_SIZE) fprintf(fp,BPOF,*(*(P+i)+j++));
				if (BTOM_MOME) fprintf(fp,BPOF,*(*(P+i)+j++));
			if (BETA) fprintf(fp,BPOF,*(*(P+i)+j++));
				if (BETA_SIZE) fprintf(fp,BPOF,*(*(P+i)+j++));
				if (BETA_BTOM) fprintf(fp,BPOF,*(*(P+i)+j++));
				if (BETA_MOME) fprintf(fp,BPOF,*(*(P+i)+j++));
			if (SICC) fprintf(fp,BPOF,*(*(P+i)+j++));
				if (SICC_SIZE) fprintf(fp,BPOF,*(*(P+i)+j++));
				if (SICC_BTOM) fprintf(fp,BPOF,*(*(P+i)+j++));
				if (SICC_MOME) fprintf(fp,BPOF,*(*(P+i)+j++));
			if (RSQR) fprintf(fp,BPOF,*(*(P+i)+j++));
			if (CLST) fprintf(fp,"%hd,",*(Q+i));
			fprintf(fp,"\n");
		}
		else {
			fprintf(fp,BPOF,**(P+i));
			fprintf(fp,"\n");
		}
	}
		
	/* close file */
	if (fclose(fp) != 0) 
		close_error(file_name);
}

/* writes a token to disk */
void tkn_write(char file_name[]) {

	char path[LLEN],line[LLEN]; /* file path and scanning string */
	
	FILE * fp; /* file pointer */
	
	/* create full file name */
	strcpy(path,DIRECTORY); 
	strcat(path,file_name);
	
	/* open file */
	if ((fp = fopen(path,"w")) == NULL) 
		open_error(file_name);
	
	/* close file */
	if (fclose(fp) != 0)
		close_error(file_name);
}