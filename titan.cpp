/* CLEAN */

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

/* R^2 Suite */
double **A_ls; /* regressor matrix for least-squares */
	double **A_temp; /* temporary 'A_ls' */
double *b_ls; /* regressand matrix for least-squares */
	double *b_temp; /* temporary 'b_ls' */
double *v_hh; /* Householder vector */
double *b_hh; /* Householder beta */
double *b_hat; /* fitted right-hand-side */

ret *R; /* returns panel */
	long R_n = 0; /* # observations */
	
	ret_mkt *R_mkt; /* market returns panel */
		short R_mkt_n = 0; /* # observations */
		
	ret_ref *R_size; /* RPRP */
	ret_ref *R_btom; 
	ret_ref *R_mome;
		short R_ref_n = 0; /* # observations */
	 
ecf *T; /* ECFP... */
	long T_n = 0; /* # observations */

	ecf *T_size; /* ...sorted by size */
	ecf *T_btom; /* ...sorted by btom */
	ecf *T_beta; /* ...sorted by beta */
	ecf *T_sicc; /* ...sorted by sicc */
	ecf *T_mome; /* ...sorted by mome */
	
	ecf *T_szbm;

ecd *D; /* event/control firm deciles */
	short D_n = 0; /* counter */

long **C; /* controls: C: MNCF x MNMC */
	short C_n = 0; /* control counter */

long *cc; /* current controls: c: MNMC x 1 ... */

double *cr; /* ...and their returns: r: MNMC x 1 */ 
	double *cr_temp; /* temporary 'cr' */

double *er; /* event returns */

node root; /* tree's root node */
	
long *B; /* bad controls */
	short B_n = 0; /* counter */
	
double **P; /* BHAR panel */
	long P_n = 0; /* # observations */
	
long *Q; /* event dates */
	long Q_n = 0; /* # observations */

long *X; /* current control glo_keys */
	long X_n = 0; 

/* gets controls */
short get_driver(ecf *);
	
/* puts BHARs */
short put_driver(ecf *);
	
/* old faithful */
int main(void) {

	long i; /* counter */
	
	/* read returns panel */
	printf("reading returns panel...\n");
	ret_read("ret.csv");
	ret_mkt_read("ret_mkt.csv");
	
	/* read event/control panel... */
	printf("reading event/control panel...\n");
	ecf_read("ecf.csv");
	
	/* ...and build it */
	printf("building event/control panels...\n");
	ecf_build();
	
	/* compute reference portfolio returns */
		
	/* compute BHARs */
	printf("computing BHARs...\n");		
		
		/* if computing an R^2, allocate memory for... */
		A_ls = f2alloc(ESTW,MNMC); /* ...the regressor matrix...*/
			A_temp = f2alloc(ESTW,MNMC); /* ...(and its temp)... */
		b_ls = f1alloc(ESTW); /* ...the regressand matrix... */
			b_temp = f1alloc(ESTW); /* ...(and its temp)... */ 
		v_hh = f1alloc(ESTW); /* ...the Householder vector... */
		b_hh = f1alloc(MNMC); /* ...the Householder beta... */
		b_hat = f1alloc(ESTW); /* ...and the fitted right-hand-side... */
		
		/* allocate memory for the bad controls */
		B = l1alloc(MNCF*MNMC);
		
		/* allocate memory for the controls */
		C = l2alloc(MNCF,MNMC);
		
			/* allocate memory for current controls */
			cc = l1alloc(MNMC);
			
			/* allocate memory for control returns */
			cr = f1alloc(MNMC);
			
				/* allocate memory for temporary 'cr' */
				cr_temp = f1alloc(MNMC);

			/* allocate memory for event returns */
			er = f1alloc(MNMC);
		
		/* allocate memory for the BHAR panel */
		P = f2alloc(MECF,MNMC+RSQR);
		
		/* allocate memory for the event dates */
		if (CLST)
			Q = l1alloc(MECF);
			
		X = l1alloc(MNMC);
		
		/* build the control tree */
		build_tree(&root,NULL,-1,-1);
		
		/* loop over candidate events */
		for (i = 0; i < T_n; i++) {
		
			/* report progress */
			if (i % 10000 == 0)
				printf("%d/%d\n", i, T_n);
				
			/* continue to next event if current event is outside of the sample window */
			if ((T+i)->year > MAXY || (T+i)->year < MINY) continue; 
		
			if ((T+i)->year>1994) continue;
		
			/* continue to next event if current event is not in desired biased sample */
			if (SML_FRMS && (T+i)->size > decl_srch('s',BDRP,(T+i)->year,(T+i)->mnth)) continue; /* small firms */
			if (LRG_FRMS && (T+i)->size < decl_srch('s',TDRP,(T+i)->year,(T+i)->mnth)) continue; /* large firms */
			if (LOW_BTOM && (T+i)->btom > decl_srch('b',BDRP,(T+i)->year,(T+i)->mnth)) continue; /* low book-to-market firms */
			if (HGH_BTOM && (T+i)->btom < decl_srch('b',TDRP,(T+i)->year,(T+i)->mnth)) continue; /* high book-to-market firms */
			if (WEA_MOME && (T+i)->mome > decl_srch('m',BDRP,(T+i)->year,(T+i)->mnth)) continue; /* weak pre-event returns */
			if (STR_MOME && (T+i)->mome < decl_srch('m',TDRP,(T+i)->year,(T+i)->mnth)) continue; /* strong pre-event returns */

			if (get_driver(T+i) == 1) continue; /* get controls */
			if (put_driver(T+i) == 1) continue; /* put controls */
		
		}

	/* output BHAR panel to disk */
	pnl_write(BPFN);
	
	/* drop the token */ 
	tkn_write("tkn.txt");
	
	/* write the ECFP to disk for analysis */
	ecf_write("ecf_out.csv");
	
	/* all ok */
	return 0;
}

/* gets controls */
short get_driver(ecf *current) {

	long twin; 

	C_n = 0; /* reset the matching char. counter */
	X_n = 0;
	/* if using [matching char.] and an insufficient number of controls are found for [matching char.], return 1 */
	
	/* size match */
	if (SIZE) 
		if ((twin = ecf_srch(T_size,current->glo_key,'n',size_pick,size_comp)) == -1)
			return 1;
		else
			*(X+X_n++) = twin;
		/* size in book-to-market decile match */
		if (SIZE_BTOM && ecf_srch(T_size,current->glo_key,'b',size_pick,size_comp) == -1) return 1;
		/* size in momentum decile match */
		if (SIZE_MOME && ecf_srch(T_size,current->glo_key,'m',size_pick,size_comp) == -1) return 1;
	/* book-to-market match */
	if (BTOM) 
		if ((twin = ecf_srch(T_btom,current->glo_key,'n',btom_pick,btom_comp)) == -1) 
			return 1;
		else
			*(X+X_n++) = twin;
		/* book-to-market in size decile match */
		if (BTOM_SIZE && ecf_srch(T_btom,current->glo_key,'s',btom_pick,btom_comp) == -1) return 1;
			/* LBT size/book-to-market match */
			if (LBT_BTOM_SIZE && ecf_srch(T_btom,current->glo_key,'l',btom_pick,btom_comp) == -1) return 1;
		/* book-to-market in momentum decile match */
		if (BTOM_MOME && ecf_srch(T_btom,current->glo_key,'m',btom_pick,btom_comp) == -1) return 1;	
	/* beta match */
	if (BETA && ecf_srch(T_beta,current->glo_key,'n',beta_pick,beta_comp) == -1) return 1;
		/* beta in size decile match */
		if (BETA_SIZE && ecf_srch(T_beta,current->glo_key,'s',beta_pick,beta_comp) == -1) return 1;
		/* beta in book-to-market decile match */
		if (BETA_BTOM && ecf_srch(T_beta,current->glo_key,'b',beta_pick,beta_comp) == -1) return 1;
		/* beta in momentum decile match */
		if (BETA_MOME && ecf_srch(T_beta,current->glo_key,'m',beta_pick,beta_comp) == -1) return 1;
		/* sic code match */
	if (SICC && ecf_srch(T_sicc,current->glo_key,'n',sicc_pick,sicc_comp) == -1) return 1;
		/* sic code in size decile match */
		if (SICC_SIZE && ecf_srch(T_sicc,current->glo_key,'s',sicc_pick,sicc_comp) == -1) return 1;
		/* sic code in book-to-market decile match */
		if (SICC_BTOM && ecf_srch(T_sicc,current->glo_key,'b',sicc_pick,sicc_comp) == -1) return 1;
		/* sic code in momentum decile match */
		if (SICC_MOME && ecf_srch(T_sicc,current->glo_key,'m',sicc_pick,sicc_comp) == -1) return 1;
	
	return 0; /* otherwise, return 0 */
}

/* puts BHARs */
short put_driver(ecf *current) {
	
	short i, j, k, bad_cgrp; /* counters and bad control group flag */
	short year = R[current->ret_key].year, mnth = R[current->ret_key].mnth; /* year and month */
	short year_temp = year, mnth_temp = mnth; /* temporary year and month */
	double er_temp, rsqr; /* temporary 'er' and R^2 */
	
	B_n = 0; /* reset the bad controls counter */ 
				
	/* reset the control and event returns */
	for (i = 0; i < MNMC; i++)
		*(cr+i) = *(er+i) = 1.;
	
	/* reset the control tree */
	reset_tree(&root,-1,-1); 
	
	/* calculate the control group ranking variables */
	if (calc_cgrv(&root,current->ret_key,year,mnth) == 1) 
		return 1;
				
	/* load the first group */
	if (next_cgrp(cc) == 1) 
		return 1;
				
	/* loop over the months in the holding period */
	for (i = 0; i < PEVW; ) {
	
		/* update temporary year and month */
		if ((mnth_temp = mnth%12+1) == 1)
			year_temp = year+1;
			
		if (MNCF > 1) {
		
			/* load returns for each control in the group */ 
			for (j = 0; j < MNMC && (bad_cgrp = ret_srch(cr_temp+j,*(cc+j)+i+1,year_temp,mnth_temp,(R+*(cc+j))->firm)) == 0; j++)
				;
			
			if (bad_cgrp) { /* if one of the controls is bad... */
				*(B+B_n++) = *(cc+j); /* ...tag it... */
				if (next_cgrp(cc) == 1) /* ...and load the next group */
					return 1; 
			}
			else { /* otherwise... */
			
				/* ...update the year and month */
				year = year_temp; 
				mnth = mnth_temp;
		
				/*
				if (ret_srch(&er_temp,current->ret_key+i+1,year,mnth,current->firm) == 0)
					for (j = 0; j < MNMC; j++)
						*(er+j) *= (1.+er_temp);
				else if (ret_mkt_srch(&er_temp,year,mnth) == 0)
					for (j = 0; j < MNMC; j++)
						return 1; //*(er+j) *= (1.+er_temp);
				else
					return 1;
				*/
		
				/* ...and update the event and controls' returns */
				for (j = 0; j < MNMC; j++) {
				
					/* update the controls' returns */
					*(cr+j) *= (1.+*(cr_temp+j));
				
					/* update the event's return */
					*(er+j) *= (ret_srch(&er_temp,current->ret_key+i+1,year,mnth,current->firm) == 0) ? (1.+er_temp) : (1.+*(cr_temp+j));

				}
				
				i++; /* update the counter */
			}
		}
		else {
		
			/* ...update the year and month */
			year = year_temp; 
			mnth = mnth_temp;

			for (j = 0; j < MNMC; j++) {
			
				if (ret_srch(cr_temp+j,*(cc+j)+i+1,year_temp,mnth_temp,(R+*(cc+j))->firm) == 0)
					*(cr+j) *= (1.+*(cr_temp+j));
				else if (ret_mkt_srch(cr_temp+j,year,mnth) == 0)
					return 1; //*(cr+j) *= (1.+*(cr_temp+j));
				else
					return 1;
			
				if (ret_srch(&er_temp,current->ret_key+i+1,year,mnth,current->firm) == 0)
					*(er+j) *= (1.+er_temp);
				else if (ret_mkt_srch(&er_temp,year,mnth) == 0)
					return 1; //*(er+j) *= (1.+er_temp);
				else
					return 1;
			}
						
			i++;
		}
	}
	
	i = 0; /* matching characteristics counter */
	j = 0;
	
	/* compute BHARs and load them into the BHAR panel */
	
		/* size match */
		if (SIZE) {
			*(*(P+P_n)+i++) = *(er++)-*(cr++);
			(T_size+*(X+j++))->used = 1;
			current->size_twin = (T_size+*(X+j-1))->glo_key;
		}
			/* size in book-to-market decile match */
			if (SIZE_BTOM) *(*(P+P_n)+i++) = *(er++)-*(cr++);
			/* size in momentum decile match */
			if (SIZE_MOME) *(*(P+P_n)+i++) = *(er++)-*(cr++);
		/* book-to-market match */
		if (BTOM) {
			*(*(P+P_n)+i++) = *(er++)-*(cr++);
			(T_btom+*(X+j++))->used = 1;
			current->btom_twin = (T_btom+*(X+j-1))->glo_key;
		}
			/* book-to-market in size decile match */
			if (BTOM_SIZE) *(*(P+P_n)+i++) = *(er++)-*(cr++);
			/* book-to-market in momentum decile match */
			if (BTOM_MOME) *(*(P+P_n)+i++) = *(er++)-*(cr++);
		/* beta match */
		if (BETA) {
			*(*(P+P_n)+i++) = *(er++)-*(cr++);
			(T_beta+*(X+j++))->used = 1;
			current->beta_twin = (T_beta+*(X+j-1))->glo_key;
		}
			/* beta in size decile match */
			if (BETA_SIZE) *(*(P+P_n)+i++) = *(er++)-*(cr++);
			/* beta in book-to-market decile match */
			if (BETA_BTOM) *(*(P+P_n)+i++) = *(er++)-*(cr++);
			/* beta in momentum decile match */
			if (BETA_MOME) *(*(P+P_n)+i++) = *(er++)-*(cr++);
			/* sic code match */
		if (SICC) *(*(P+P_n)+i++) = *(er++)-*(cr++);
			/* sic code in size decile match */
			if (SICC_SIZE) *(*(P+P_n)+i++) = *(er++)-*(cr++);
			/* sic code in book-to-market decile match */
			if (SICC_BTOM) *(*(P+P_n)+i++) = *(er++)-*(cr++);
			/* sic code in momentum decile match */
			if (SICC_MOME) *(*(P+P_n)+i++) = *(er++)-*(cr++);
		
		/* compute R^2s and load them into the BHAR panel */
		if (RSQR && mrsq(*(P+P_n)+i++,cc,current->ret_key,(R+current->ret_key)->year,(R+current->ret_key)->mnth) == 1)
			return 1;
			
		/* load year and month into the BHAR panel */
		if (CLST)
			*(Q+Q_n++) = 100*current->year+current->mnth;
		
		P_n++; /* update the BHAR panel counter */
		cr -= MNMC; er -= MNMC; /* reset the pointers */
	
	/* all ok */
	return 0;
}