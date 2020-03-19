/* comparison functions */

#include <stdlib.h>
#include <math.h>

#include "glib.h"

int ymp_comp(const void *p1, const void *p2) {

	const ymp *year_mnth1 = (const ymp *) p1;
	const ymp *year_mnth2 = (const ymp *) p2;
	
	if (year_mnth1->year < year_mnth2->year) 
		return -1;
	else if (year_mnth1->year == year_mnth2->year)
		if (year_mnth1->mnth < year_mnth2->mnth)
			return -1;
		else if (year_mnth1->mnth == year_mnth2->mnth)
			return 0;
		else
			return 1;
	else
		return 1;
}

int rkey_comp1(const void *p1, const void *p2) {

	const int *a1 = &(((const ret *) p1)->rkey);
	const int *a2 = &(((const ret *) p2)->rkey);

	return (*a1>*a2)-(*a1<*a2);
}

int rkey_comp2(const void *p1, const void *p2) {

	const int *a1 = &(((const ecf *) p1)->rkey);
	const int *a2 = &(((const ecf *) p2)->rkey);

	return (*a1>*a2)-(*a1<*a2);
}