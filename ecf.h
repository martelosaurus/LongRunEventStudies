#ifndef TABLES_H

	#define TABLES_H

	typedef struct ecf_struct {
		
		/* keys */
		long ret_key; /* key for returns */
 		long glo_key; /* key for ECFP sorted by firm */
		long loc_key; /* key for sorted ECFP */
		
		/* ids */
		short year;
		short mnth;
		long firm;
		short excd;
		
		/* factors */
		double size;
		double btom;
		double beta;
		short sicc;
		double mome;
		
		/* decile ranks */
		short size_rank;
		short btom_rank;
		short mome_rank;
		
		/* twins */
		short used;
		long size_twin;
		long btom_twin;
		long beta_twin;
		
	} ecf;

	typedef struct ecd_struct {
	
		/* ids */
		short year;
		short mnth;
		short rank;
	
		/* deciles */
		double size;
		double btom;
		double mome;
		
	} ecd;

	/* BUILD FUNCTIONS */
	
		/* creates sorted copies of ECFP and computes decile ranks */
		void ecf_build(void);
		
		/* copies contents of one ECFP to another */
		long ecf_copy(ecf *, ecf *, short, short, short);
		
		/* creates sorted copies of ECFP */
		void ecf_sort(ecf *, int (*comp)(const void *, const void *));
	
	/* DECILE FUNCTIONS */
		
		/* computes deciles */
		void ecf_decls(ecf *, short, void (* pick)(ecd *, const ecf *));
		
		/* computes decile ranks for each observation in ECFP */
		void ecf_ranks(ecf *);
	
		/* searches decile object for desired decile */
		double decl_srch(char, short, short, short);	
	
	/* SEARCH FUNCTIONS */
	
		/* searches the specified ECFP for a set of controls for the specified event */ 
		long ecf_srch(ecf *, long, char df, void (* pick)(ecf *, const ecf *), int (* comp)(const void *, const void *));
		
		/* subroutine of ecf_srch: performs local-linear search */
		long ecf_loc_srch(ecf *, long, long, char df);
	
#endif