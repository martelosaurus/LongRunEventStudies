#include <stdlib.h>

#include "levlib.h"

grp *stack_pop(stack *S) {
    if (S->n == 0)
        return NULL;
    else
        return S->cgrp+(--S->n); 
}

void stack_sort(stack *S) {
    qsort(S->cgrp,S->N,sizeof(grp),cgrp_cgrv_comp);
}

int cgrp_cgrv_comp(const void *p1, const void *p2) {
    
    const double *a1 = &(((const grp *) p1)->rsqr);
	const double *a2 = &(((const grp *) p2)->rsqr);

    /* lowest to highest */
	return (*a1>*a2)-(*a1<*a2);
}