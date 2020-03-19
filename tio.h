#ifndef TIO_H

	#define TIO_H
	
	/* reads the ECFP from disk */
	void ecf_read(char file_name[]);
	
	/* writes the ECFP to disk */
	void ecf_write(char file_name[]);
	
	/* reads the returns panel from disk */
	void ret_read(char file_name[]);
	
	/* reads the market returns panel from disk */
	void ret_mkt_read(char file_name[]);
	
	/* writes the BHAR panel to disk */
	void pnl_write(char file_name[]);
	
	/* writes a token to disk */
	void tkn_write(char file_name[]);
	
#endif 