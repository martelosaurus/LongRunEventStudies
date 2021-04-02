#ifndef TPS_H

	#define TPS_H
	
	/* Lexicon 
	
		ECFP: event/control firm panel: the base panel is the one sorted by firm
		BHAR: buy-and-hold abnormal returns
	
		size: market capitalization
		btom: book-to-market ratio
		beta: market model beta
		sicc: 2-digit SIC code
		mome: pre-event holding return
		
		MECF: maximum number of observations from ecf.csv
		MRET: ----------------------------------- ret.csv

		MNCF: maximum number of backup control firms plus one
		BHOR: backward horizon (in months) for computing beta and pre-event returns
		PEVW: forward horizon (in months) for computing buy-and-hold abnormal returns
		BDRP: bottom decile rank position: values < D[BDRP]. are in the lowest decile
		TDRP: top decile rank position: values > D[TDRP]. are in the highest decile
		LLEN: maximum line length
		BPOF: BHAR panel output format
		MNMC: maximum number of matching characteristics
		MINY: minimum year
		MAXY: maximum year
		MNOM: maximum number of months
		MNRP: maximum number of reference portfolios
		
		MACH: matching characteristics
		MRSQ: maximal R^2
		LBTC: Lyon, Barber and Tsai (LBT) control firms 
		LBTP: LBT reference portfolios
		
		CGRV: control group ranking variable: either a sorting rank or an R^2
		
	*/
	
	/* Notes:
	
		- All returns are of the form Return = (Income + Price_{t+1} - Price_{t}) / Price_{t}.
		- A function return value of 0 indicates success while a return value of 1 indicates failure.
	*/
	
	/* method */
	#define MACH 1
	#define MRSQ 0
	#define LBTC 0
	#define LBTP 0
	
	/* controls and matching chars. */
	#define MNCF 1
	#define MNMC 3
	#define MNRP 50
	
	/* windows and dates */
	#define ESTW 60
	#define MOMW 6
	#define PEVW 12
	#define MINY 1963
	#define MAXY 2014
	#define MNOM 12*(MAXY-MINY+1)
	
	/* reporting */
	#define RSQR 0
	#define CLST 0
	
	/* sample */
	#define SML_FRMS 0
	#define LRG_FRMS 0
	#define LOW_BTOM 0
	#define HGH_BTOM 0
	#define WEA_MOME 0
	#define STR_MOME 0
		
	/* matching */
	#define SIZE 1
		#define SIZE_BTOM 0
		#define SIZE_MOME 0
	#define BTOM 1
		#define BTOM_SIZE 0
			#define LBT_BTOM_SIZE 0 
		#define BTOM_MOME 0
	#define BETA 1
		#define BETA_SIZE 0
		#define BETA_BTOM 0
		#define BETA_MOME 0
	#define SICC 0
		#define SICC_SIZE 0
		#define SICC_BTOM 0
		#define SICC_MOME 0
	#define MOME 0
		#define MOME_SIZE 0
		#define MOME_BTOM 0
		
	/* paths */
	#define DIRECTORY "E:\\Users\\joma1199\\LRES\\"  
	#define BPFN "exp_pnl.csv"
		
	/* observations */
	#define MECF 1799697
	#define MRET 2895953
	
	/* decile ranks */
	#define BDRP 0
	#define TDRP 8
	
	/* output */
	#define LLEN 1000
	#define BPOF "%3.12f,"
	
#endif