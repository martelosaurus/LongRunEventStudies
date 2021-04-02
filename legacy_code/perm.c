#include <stdlib.h>
#include <stdio.h>
#include "perm.h"

void perm_get(perm_list *X, perm_node *crnt) {
    
	int k, *path = X->L[X->n];
	if (crnt->next == NULL) {
		for (; crnt->last != NULL; crnt = crnt->last)
		    *(path++) = crnt->tag;
		X->n++;
	}
    else
	    for (k = 0; k < X->k; k++)
		    perm_get(X, crnt->next+k);
}

void perm_put(perm_list *X, perm_node *crnt, perm_node *last, int m, int k) {
		 
    crnt->tag = k;
	crnt->last = last;
    if (m == X->m)
		crnt->next = NULL;
    else {
		crnt->next = (perm_node *) malloc(X->k*sizeof(perm_node));
		for (k = 0; k < X->k; k++)
		    perm_put(X, crnt->next+k, crnt, m+1, k);
	}
}

void perm_pop(perm_list *X, perm_node *crnt) {
	
	int k;
	
	printf("pop\n");

	if (crnt->next == NULL)
		free(crnt);
	else 
		for (k = 0; k < X->k; k++)
			perm_pop(X, crnt->next + k);
}

int perm_power(int b, int e) {
		
    int p;

	for (p = 1; e > 0; --e)
	    p *= b;
    return p;
}

int **perm(int k_max, int m_max) {

	int i, n, **r;
	perm_node root;
	perm_list X;
	
	n = perm_power(k_max, m_max);

	X.k = k_max;
	X.m = m_max;
	X.n = 0;
	X.L = (int **) malloc(n*sizeof(int *));
	
	for (i = 0; i < n; i++)
		X.L[i] = (int *) malloc(X.m*sizeof(int));
	
	perm_put(&X, &root, NULL, 0, 0);
	perm_get(&X, &root);
	/* perm_pop(&X, &root); */
	
	r = X.L;
	
	return r;
}