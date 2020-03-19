#ifndef RET_H

	#define RET_H
	
	typedef struct ret_struct {
	
		// keys
		long ret_key;
		
		// return id
		long firm;
		short year;
		short mnth;
		
		// return
		double fret;
		
	} ret;
	
	typedef struct ret_mkt_struct {
		
		// return id
		short year;
		short mnth;
		
		// return
		double mret;
	
	} ret_mkt;
	
	typedef struct ret_ref_struct {
	
		// return id
		short rank;
		short year;
		short mnth;
		
		// return
		double pret;
	
	} ret_ref;
	
	/* computes an event/control's beta */
	short calc_beta(double *, long);
	
	/* computes an event/control's momentum */
	short calc_mome(double *, long);
	
	/* computes the reference portfolio returns */
	void ecf_build(void);
	void ref_calc(ret_ref *);
	
	/* searches the returns panel for the specified event/control and updates firm/market returns */ 
	short ret_srch(double *, long, short, short, long);
	
	short ret_mkt_srch(double *, short, short);
	
	short ret_ref_srch(char, double *, short, short, short, short);
	
#endif