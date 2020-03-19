/* titan trees and climbers */

#include <stdio.h>
#include <stdlib.h>

#include "ecf.h" 
#include "ret.h" 
#include "tio.h"
#include "tls.h" 
#include "tps.h" 
#include "ttc.h"
#include "tuf.h"

extern double **A_ls; /* regressor matrix for least-squares */
	extern double **A_temp; /* temporary 'A_ls' */
extern double *b_ls; /* regressand matrix for least-squares */
	extern double *b_temp; /* temporary 'b_ls' */

extern ret *R; /* returns panel */
	extern long R_n; /* # observations */

extern long *B; /* bad controls */
	extern short B_n; /* counter */
	
extern long **C; /* controls: C: MNCF x MNMC */
	extern short C_n; /* control counter */
	
extern node root; /* tree's root node */

/* build the tree */ 
/* for the root, set j = -1 */
void build_tree(node *current, node *last, short j, short i) {

	short k; /* counter */

	current->rank = i; /* set the control rank */
	current->last = last; /* a pointer to the parent node */
	
	if (j == MNMC-1) { /* if at the top of the tree... */
		
		current->next = NULL; /* ...indicate so */
		
		/* allocate memory for the paths */
		current->key_path = l1alloc(MNMC);
		current->rnk_path = s1alloc(MNMC);
		
		/* build the rank path */
		build_rnk_path(current,current->rnk_path);
	}
	else { /* otherwise... */
		
		/* allocate memory for the child nodes... */
		if ((current->next = (node *) malloc(MNCF*sizeof(node))) == NULL)
			general_error("error: mem_all: build_tree");
		
		/* ...and build them */
		for (k = 0; k < MNCF; k++)
			build_tree(current->next+k,current,j+1,k);
	}
}

/* reset the tree (for a new event) */
/* for the root, set j = -1 */
void reset_tree(node *current, short j, short i) {
	
	short k; /* counter */
	
	current->ctrl = (j != -1) ? *(*(C+i)+j) : -1; /* set the control associated with the current node */

	if (current->next == NULL) /* if at the top of the tree... */
		build_key_path(current,current->key_path); /* ...build the key path */
	else /* otherwise... */
		for (k = 0; k < MNCF; k++) /* ...continue toward the root */
			reset_tree(current->next+k,j+1,k);
}

/* destroy the tree (inactive) */
void destroy_tree(node *current) {
	short i;
	if (current->next == NULL) {
		free(current->key_path); 
		free(current->rnk_path);
	}
	else {
		for (i = 0; i < MNCF; i++)
			destroy_tree(current->next+i);
		free(current->next);
	}
}

/* build the key path from the current node to the root */
void build_key_path(node *current, long *key_path) {

	if (current->last != NULL) { /* if the node is not the root... */
		*key_path = current->ctrl; /* ...add the node's ctrl to the key_path... */
		build_key_path(current->last,++key_path); /* ...and continue toward the root */
	}
}

/* build the rank path from the current node to the root */
void build_rnk_path(node *current, short *rnk_path) {

	if (current->last != NULL) { /* if the node is not the root... */
		*rnk_path = current->rank; /* ...add the node's rank to the rnk_path... */
		build_rnk_path(current->last,++rnk_path); /* ...and continue toward the root */
	}
}

/* calculate the control group's ranking variable */
short calc_cgrv(node *current, long evnt, short year, short mnth) {
	
	/* evnt: RET_KEY of event */
	
	short i; /* counters */
	
	if (current->next == NULL) { /* if at the top of the tree... */
		if (MACH) /* ...and using matching characteristics... */
			current->cgrv = mach(current->rnk_path);  /* ...return matching characteristics CGRV */
		else /* ...and using maximal R^2... */
			if (mrsq(&current->cgrv,current->key_path,evnt,year,mnth) == 1) /* ...return matching characteristics CGRV */
				return 1; 
	}
	else /* otherwise, continue toward the top of the tree */
		for (i = 0; i < MNCF; i++) 
			if (calc_cgrv(current->next+i,evnt,year,mnth) == 1)
				return 1;
			
	return 0; /* all ok */
}

/* wrapper for best_cgrv */
short next_cgrp(long *key_path) {
	return (best_cgrv(&root,key_path) == NULL) ? 1 : 0; /* run the driver */
}

/* searches for next best control group */
double *best_cgrv(node *current, long *key_path) {
	
	int i, flag = 0; /* counter and flag */
	double *temp_val, best_val = -1.; /* temporary CGRV and best CGRV */
	
	if (current->next == NULL) { /* if at top of the tree... */
		/* return a NULL pointer if the "key_path" (control group) is bad and a pointer to the CGRV otherwise */
		return (is_bad(current->key_path)) ? NULL : &current->cgrv; 
	}
	else { /* otherwise... */
		for (i = 0; i < MNCF; i++) { /* ...find the maximum CGRV from the current node's children */
			temp_val = best_cgrv(current->next+i,key_path+1); /* a pointer to the child's best CGRV (or NULL) */
			if (temp_val != NULL && *temp_val <= best_val) { /* if it's not NULL and at least as big as the best CGRV... */
				best_val = *temp_val; /* ...update the best CGRV... */
				*key_path = (current->next+i)->ctrl; /* ...update the "key_path" (control group)... */
				flag = 1; /* if a non-NULL CGRV is found, set the flag to one */
			}
		}
		/* if a non-NULL CGRV was found, return a pointer to the best CGRV; otherwise, return NULL */
		return (flag) ? &best_val : NULL; 
	}
}

/* checks whether a "key_path" (control group) is bad */
short is_bad(long *key_path) {

	int i, j; /* counters */
	
	for (i = 0; i < B_n; i++) /* loop over the bad controls */
		for (j = 0; j < MNCF; j++) /* loop over the controls in the group */
			if (*(key_path+j) == *(B+i)) /* if there's a match... */
				return 1; /* ...return an error code */
				
	return 0; /* all ok */
}

/* calculates matching characteristics ranking variable */
double mach(short *rnk_path) {
	
	short i; /* counter */
	double cgrv = 0.; /* initialize the CGRV */
	
	/* the CGRV as the sum of the inverse ranks of all controls in the group */
	for (i = 0; i < MNCF; i++)
		cgrv += MNCF-*(rnk_path++);
		
	return cgrv; /* return the CGRV */
}

/* calculates maximal R^2 ranking variable */
short mrsq(double *rsqr, long *key_path, long evnt, short year, short mnth) {

	/* evnt: RET_KEY of event */

	short i, j, year_temp, mnth_temp; /* counters and temporary year and month */
	static long prev_evnt = -1; /* previous event */
	static short prev_ctrl_alloc = 0; /* alloc flag */
	static long *prev_ctrl; /* previous controls */
	
	if (!prev_ctrl_alloc) { /* if the previous control array has not been alloced */
		if ((prev_ctrl = (long *) malloc(MNMC*sizeof(long))) == NULL) /* ...allocate memory for it... */
			general_error("error: mem_all: mrsq");
		for (i = 0; i < MNMC; i++) /* ...initialize each entry to -1... */
			*(prev_ctrl+i) = -1;
		prev_ctrl_alloc = 1; /* ...and reset the flag */
	}
	
	/* for each control in the group, copy it's returns into the regressand matrix */
	for (j = 0; j < MNMC; j++) { 
	
		if (*(key_path+j) == *(prev_ctrl+j)) /* if the control was previously used in the same position, continue */
			continue;
	
		/* initialize and loop over the months in the sample window */
		year_temp = year; mnth_temp = mnth; 
		for (i = 0; i < ESTW; i++) {
		
			/* update year and mnth */
			if ((mnth_temp = 12-(13-mnth_temp)%12) == 12)
				year_temp--;
		
			/* copy the control's returns into the regressand matrix */
			if (ret_srch(*(A_temp+i)+j,*(key_path+j)-i-1,year_temp,mnth_temp,(R+*(key_path+j))->firm) == 1)
				return 1;
		}
		
		*(prev_ctrl+j) = *(key_path+j); /* update the previous controls */
	}
	
	
	if (evnt != prev_evnt) { /* if the current event is different than the previous event... */
	
		/* initialize and loop over the months in the sample window */
		year_temp = year; mnth_temp = mnth;
		for (i = 0; i < ESTW; i++) {
			
			/* update year and mnth */
			if ((mnth_temp = 12-(13-mnth_temp)%12) == 12)
				year_temp--;
		
			/* copy the event's returns into the regressor matrix */
			if (ret_srch(b_temp+i,evnt-i-1,year_temp,mnth_temp,(R+evnt)->firm) == 1)  
				return 1;
		}
		
		prev_evnt = evnt; /* update the previous event */
	}
	
	*rsqr = calc_vrnc(evnt); /* compute R^2 */
		
	return 0; /* all ok */
}