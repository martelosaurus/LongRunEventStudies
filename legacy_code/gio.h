#ifndef KIO_H

	#define KIO_H
	
	/* missing beta code */
	#define MBC -99.
	#define LEN 100
	#define DIRECTORY "\\\\BusApollo\\Users\\joma1199\\grendel\\"
	
	int get_ret(char *);
	int get_ecf(char *);
	int put_ecf(char *);
	int put_pnl(char *, int);
	FILE *open_file(char *, char *);
	
#endif