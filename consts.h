#define FIRLEN 4	//Hard-coded for each filter. TODO - make this more flexible, e.g. read in a config file, or pass from a C shell, or include a header
					//in the filter dat file
#define MAXFILTERLENGTH 64  //Assume no filter will be longer than this. CAVEAT!
#define TAPLENGTH (MAXFILTERLENGTH + 2) //Even if filter has max length, there is room for 2 new samples to be written without disturbing ongoing convolution
#define DELAY_SIZE  TAPLENGTH //#def BASE_TAPLENGTH in fir_coeff.h, =8, Now making d[] always big enough (pseudo-)I1 gets repositioned to three samples behind I0 before each
		//convolution

//flag bits
#define IN_READY 0	//low bit in sam (e.g.) signals when new samples are in
#define ONE_DONE 1	//ISR checks this to see if R7 contains one 16-bit sample already.

#define DONE_LP 2	//set on doing LP (scaling, approximation) filter, thereupon do HP (wavelet, detail) filter
#define BYE 15 //If set main will exit 

#define NLEVELS 10  //How many resolution levels, = no.decimations = no.of delay lines to allocate

#define SAMPSIZE 2	//16-bit samples