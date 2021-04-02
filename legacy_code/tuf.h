#ifndef TUF_H
	
	#define TUF_H
	
	/* computes the sample covariance of x and y (see definition for additional notes) */
	double cov(double *x, double *y); 
	
	/* reports... */
	void general_error(char *error); /* ...general error */
	void open_error(char *file_name); /* ...error in opening a file */
	void close_error(char *file_name); /* ...error in closing a file */
	
	/* allocates memory for... */
	short *s1alloc(long); /* ...an m x 1 array of shorts */
	long *l1alloc(long); /* ...an m x 1 array of longs */
	double *f1alloc(long); /* ...an m x 1 array of floats */
	long **l2alloc(long, long); /* ...an m x n array of longs */
	double **f2alloc(long, long); /* ...and m x n array of floats */
	
	/* assigns factor value from 'from' to 'to' */
	void size_decl_pick(ecd *to, const ecf *from);
	void btom_decl_pick(ecd *to, const ecf *from);
	void mome_decl_pick(ecd *to, const ecf *from);
	
	/* assigns factor value from 'from' to 'to' */
	short rank_pick(ecf *x, char mach);
	void size_pick(ecf *to, const ecf *from);
	void btom_pick(ecf *to, const ecf *from);	
	void beta_pick(ecf *to, const ecf *from);
	void sicc_pick(ecf *to, const ecf *from);
	void mome_pick(ecf *to, const ecf *from);
	
	/* firm comparison function... */ 
	int ecf_ret_key_comp(const void *p1, const void *p2); /* ...for returns panel */
	int ret_ret_key_comp(const void *p1, const void *p2); /* ...for ECFP */
	
	/* factor comparison functions */
	int size_comp(const void *p1, const void *p2);
	int btom_comp(const void *p1, const void *p2);
	int beta_comp(const void *p1, const void *p2);
	int sicc_comp(const void *p1, const void *p2);
	int mome_comp(const void *p1, const void *p2);
	int szbm_comp(const void *p1, const void *p2);
	
#endif
