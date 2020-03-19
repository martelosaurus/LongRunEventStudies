#ifndef KLIB_H

	#define KLIB_H
	
	#define RET 2895953
	#define ECF 5000
	#define MNCF 3
	#define MNMC 3
	#define MLSR 10
	#define FRWD 12
	#define BKWD 60
	#define MAXY 1994
	#define MINY 1963
	#define NYSE 1
	#define MNOM 12*(MAXY-MINY+1)
	#define QNTL 10
	#define SIZE 0
	#define ECF_OUT 1
	
	/* returns */
	typedef struct ret_struct {
		
		/* keys */
		int rkey;
	
		/* ids */
		int firm;
		int year;
		int mnth;
		
		/* returns */
		double rtrn;
		double mrkt;
		
	} ret;
	
	/* event/control firm */
	typedef struct ecf_struct {
	
		/* keys */
		int fkey; 
		int rkey; 
		
		/* ids */
		int firm; 
		int excd;
		int year; 
		int mnth;
		
		/* factors and ranks */
		double *fctr;
		int *rank;
		
		/* twins */
		double *twin1;
		double *twin2;
	
	} ecf;
	
	/* hash table */
	typedef struct htab {

		int a;
		int b;

	} htab;

	/* year-month-pointer */
	typedef struct ymp_struct {
	
		/* ids */
		int year;
		int mnth;
		
		/* pointers */
		ecf *pntr;
		htab ***hash;
		
		/* counters */
		int n;
		
	} ymp;
	
	/* matching factor specs */
	typedef struct mfs_struct {
		
		/* lower*fctr2 <= fctr1 <= upper*fctr2 */
		int fctr1;
		int fctr2;
		double fctr2_upper;
		double fctr2_lower;
		
		/* quantile matching factor */ 
		int fctr3;
		int fctr3_upper;
		int fctr3_lower;
	
	} mfs;
	
	/* biased sample specs */
	typedef struct bss_struct {
		
		int fctr;
		int fctr_upper;
		int fctr_lower;
		
	} bss;
	
	/* event/control firm group */ 
	typedef struct group_struct {
		
		ecf **firm; /* control group stocks */
		double rsqr; /* control group r-squared */
		double *beta; /* control group betas */
		
	} grp;
	
	/* event/control firm group stack */
    typedef struct stack_struct {
        
		int n; /* current size */
		int N; /* maximum size */
        grp *cgrp;
    
	} stack;
	
	/* glib1: event/control firm functions */
	int new_ecf(void);
	int calc_beta(double *, ecf *);
	int calc_rank(ymp *, int, int);
	int ecf_comp(const void *p1, const void *p2);
	int htab_compb(const void *p1, const void *p2);
	
	/* glib2: search functions */ 
	ret *ret_srch(ecf *, int, int, int);
	int ecf_srch(ecf *[], ecf *, ymp *, mfs *);
	
	/* glib3: memory allocators */
	int alloc1(void);
	int alloc2(int);
	int alloc3(int);
	
	/* glib4: comparison functions */
	int ymp_comp(const void *p1, const void *p2);
	int rkey_comp1(const void *p1, const void *p2);
	int rkey_comp2(const void *p1, const void *p2);
	
	/* glib5: stack functions */
    grp *stack_pop(stack *);
	void stack_sort(stack *);
	int cgrp_cgrv_comp(const void *, const void *);
	
#endif