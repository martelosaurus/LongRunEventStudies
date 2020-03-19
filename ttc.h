#ifndef TTC_H

	#define TTC_H
	
	/* pointer to node structure */
	typedef struct node_struct *Node;

	/* node structure */
	typedef struct node_struct {
		long  ctrl; /* control ret_key */
		short  rank; /* control group rank */
		double  cgrv; /* control group ranking variable */
		long  *key_path; /* path of ctrls from node to root */
		short *rnk_path; /* path of ranks form node to root */
		Node last; /* pointer to the parent node*/
		Node next; /* pointer to the first child */
	} node;
	
	/* BUILD FUNCTIONS */
	
		/* build the tree */ 
		void build_tree(node *, node *, short, short);
		
		/* reset the tree (for a new event) */
		void reset_tree(node *, short, short);
		
		/* destroy the tree */
		void destroy_tree(node *);
		
		/* build paths from the current node to the root... */
		void build_key_path(node *, long *); /* ...using the ret_keys of the controls... */
		void build_rnk_path(node *, short *); /* ...and the tree ranks of the controls... */
		
		/* calculate the control group's ranking variable */
		short calc_cgrv(node *, long, short, short);
	
	/* SEARCH FUNCTIONS */
	
		/* wrapper for best_cgrv */
		short next_cgrp(long *);
		
		/* searches for next best control group */
		double *best_cgrv(node *, long *);
		
		/* checks whether a "path" (control group) is bad */
		short is_bad(long *);
	
	/* RANKING VARIABLE CALCULATORS */
	
		/* calculates matching characteristics ranking variable */
		double mach(short *);
		
		/* calculates maximal R^2 ranking variable */
		short mrsq(double *, long *, long, short, short);
	
#endif