

/************
Now using delay line of Mmax + 2 samples, so I1 ends up where I0 is, before next 2 samples are written to [I0++]. Need to back I1 by 1 sample.
Here's an illustration for the case of Mmax = 8 (in practice much larger), showing the transition between getting out samples y[8], y[9] and initiating
convolutions to get y[6], y[7]

Before convoluting to get y[8], y[9]
                                                                                ----------- used to get y[9] ---------------
                                                        			----------- used to get y[8] -----------------
                                                          
 |    x[0]     |    x[1]     |    x[2]     |    x[3]     |  x[4]   |    x[5]    |    x[6]   |   x[7]   |   x[8]   |   x[9]   | 
 ^                                                                                           ^
 I0                                                                                        "I1"
                                
"I1" is a pseudonym for a cached value of I0 now used as the moving pointer in the convolution where I previously used I1, making another I pointer
available elsewhere.

After calculating y[8], y[9], and at some point during this time loading x[10], x[11] at [I0++]:

 |    x[10]     |    x[11]     |    x[2]     |    x[3]     |  x[4]   |    x[5]    |    x[6]   |   x[7]   |   x[8]   |   x[9]   | 
                                ^                                     ^
                                I0                                   "I1"
                                
  The last operation using I1 was R0.H = W[I1--], getting x[4]

*****************/

/****
There are two outputs from each level, detail and approximation (high-pass and low-pass). The low-pass output goes
into the delay line of the next decimation level. The high-pass could go either to file if we want to store the transform
or also to a circular buffer, where a certain amount of history will be available for on-the-fly processing before being
fed to the reconstruction edacsac (reversed cascade). Let us assume the latter for the present, we are looking to enable
real-time processing of the wavelet coeffs and reconstruction.

Now at resolution level r. At r-1 two samples were written in to our current delay line at, let's say I0.

For all level we have
	- the high-pass filter coeff array hh
	- the low-pass filter coeff array
	 
At each level there are 4 circular buffers: 
0 	- the rth delay line, alias input buffer
1 	- belay that, the output from this level's lp filter is cached until the next level,then written to deli[r+1]- the r+1th delay line, alias low-pass (scaling function) output buffer
2	- the high-pass output buffer (wavelet coeffs) accessible for processing before reconstruction

To set up the rth level we must get B (fixed), L (fixed) and I (changing) values for each of the above.
Let's have these in a frame allocated with LINK (size will be e.g. <sizeof B+L+I> * 4 + 4, last 4 for 2-bit flags, or dedicated an Rreg to hold flags.
See LINK/UNLINK in help.

the two-bit flags (2 bits for each level, 0 means no new samples in, 1 means one sample in waiting for second, 2 means ready to process)
-OR, newthink: count input samples, N

_______________________________________________________________________________________________________________________________________________
|Do level 0 (where the input sample pairs come from the outside world via DMA) in the DMA ISR.    or maybe not, just write the captured
|samples to level-0 deli and update flags. RLOOP checks flags and decides when to run 0-level convolution                                |
|______________________________________________________________________________________________________________________________________________|
N.B., this means the next deli, at deli[TAPLENGTH*1], will get new samples written to it willy-nilly, on every ISR.Therefore must have separate 
circular index pointers.  in fact let us adhere to this principle throughout the cascade, viz:
	set B0 to the currently processing deli[r][0]
	get saved I0 pointing at last inserted sample pair
	set B1 to deli[r+1][0]
	get saved I1 pointing to where current level convolution LP output should go
	set B2 to waco[r]
	get saved I2 pointing to where HP output of current convolution should go
	

	
Make frame space for the two-bit "ready to process" flags, 32 bits, e.g. at [FP] itself. No, cannot BITTST/BITSET on [FP].
Maybe nevertheless dedicate an Rreg to this.
OR, use byte-length counters in the frame, make it 16 for good alignment

**More thoughts about cascade counters**
We know that level r will have two new samples, i.e. be ready for the next convolution, after 2^(r+2) samples have been processed in the input stream.

hh[], hl[] are always the same, L2 = M, and after the last convolution at level r-1, or at least before we start this level,
we will have done R2 = [I2++], so that once again R2.L contains h[0], R2.H h[1].

See "saving and resuming loops" for info on what to do if LTx, LBx, LCx are used in interrupts and in "user mode" code 

TIMING CONSIDERATIONS
All levels are processed in sequence, but the sequence (RLOOP) can be interrupted at any time by the capture ISR. At this time the ISR caches two samples (on
successive IRQ's) in R7. When two samples are ready in R7 the ISR must flag this to RLOOP. 
The RLOOP will not deal with those two new samples until it has run its course and started from the top (whereupon it will spin until it sees the flag).
This means that any lower levels pending will get done (generating detail and approximation coefficients, the later fed to the next lower level) before level
0 is dealt with. If we can't keep up and new samples are put into R7 before the old ones are dealt with then we get rubbish, but being fast enough to keep up
is a sine qua non of the whole algorithm.
Instead of a flag we could just keep a count of samples captured, let the counter be sampcnt. If !(sampcnt&1) then two samples have been captured, 
But for now, after two samps are captured a flag bit is set in <flags>. When RLOOP gets to check this it disables interrupts, copies R7 to deli[0][] and runs
level 0 filters on deli[0][].


****/

#include "consts.h"

.section dlines
.align 4
.byte2 deli[TAPLENGTH*NLEVELS]; //Allocate space for NLEVELS delay lines length TAPLENGTH, two-byte samples
.byte2 waco[TAPLENGTH*NLEVELS]; //where we keep wavelet coeffs (hi-pass output), for a while, so we can may do stuff with them, then reconstruct

//May be silly, since at least deli's must be at least as long as filters, REVISIT

.section firs
.VAR hh[MAXFILTERLENGTH] = "db2hpf.dat"; //Allocate as for largest allowed filter, setting L2 will make it circular at the correct length
.VAR hl[MAXFILTERLENGTH] = "db2lpf.dat"; //N.B. Help is wrong, these files must be comma-delimited, else craziness (tells me db2hpf assigns only 3 values)

.section misc
.byte4 flags = 0; //We'll never need 32 bit-flags, using byte4 for alignment


//....................................................

.align 4 //Had .align 8 here but then some addressing by multiples of 4 went cockeyed
.section L1_code;

_main:

/*** allow frame space for "structs" which look like
	{
		[totalsamps						4 bytes, count of total number of samples captured, only at level 0]
		insampcnt						4 bytes, counts how many samps fed in from level above waiting to be processed, should never exceed 2
		drr0 - address of deli[r]; 		4 bytes, origin of circular buffer we are convoluting 
		i0r - current I0[r];			4 bytes, where to write next two samples/sliding pointer in deli[r]
		wcr0 - address of waco[r];		4 bytes, origin of circular buffer for high-pass output, if we want to process on the fly
		i3r - current I3[r];			4 bytes, sliding pointer in waco[r], where we save high-pass output samples/wavelet coeffs
	}
	
	Call that rstruct (just a name, structure for resolution level r, typ)
		
or perhaps, to avoid getting deli[r+1], cache output from r in a reg or on stack, flag data ready somehow, write in next level when its turn comes.

Viewing what we have in the frame as an array of rstructs, in C speak

	rstruct rstructs[NLEVELS]
		
with rstructs = SP after the frame is constructed, so that SP is then pointing at rstructs[0]
	
ARITHMETIC OPERATIONS INSTRUCTIONS to see what you can (hence also can't) do

**********/
#define RFRAME_SIZE_B 20 //Create a "struct" this size in the frame for each level, for saving delay line sample count, addresses and pointers therein
#define HOP_OVER (RFRAME_SIZE_B - 4) //Hop over this much to get to struct for next frame if don't have two new samples to process
#define RLEV 10	//Number of resolution levels. Assuming Nyquist ~22 kHz, lowest res level will correspond to ~ 20 Hz

//	LINK NLEVELS*RFRAME_SIZE_B + 4; //assuming the six elements/addresses above for each level - no, create the frame push by push
//After this FP points at the frame space top address, SP at the bottom
//Can just as well make FP=SP, then push all required items on the stack, or even don't bother with FP, just pop items off the stack as needed. Of course in
//that case SP needs to saved by any ISR which might alter it, but ISR's use the supervisor stack, so for now no worries
	
//no longer needed, P5 = SP; //Start from bottom of frame, as it seems we can't subtract from pointers ????
//instead save current SP as FP just in case
	FP = SP;

	L0 = TAPLENGTH<<1; //x2 as 16-bit samples
	M0 = SAMPSIZE; //Use this to get from I0 pointing at address after the last two samples added to a deli, call then x[n], x[n+1], to the address
			//of x[n+1]
			
//Set up circular buff 2 for the low-pass filter coeffs
	L2 = FIRLEN<<1;	//Number of filter coeffs, x2 as 2 bytes long, possibly extended to even, makes h[] circular, of the correct length, <= MAXFILTERLENGTH
	R2.L = hh;
	R2.H = hh;
	I2 = R2;
	B2 = R2;
	
//Set up LBI 3 for high-pass
	L3 = FIRLEN;
	R3.L = hl;
	R3.H = hl;
	I3 = R3;
	B3 = R3;	 
	
//Set up address structs in frame
	P5 = TAPLENGTH<<1;
	P2=NLEVELS;
	P4.L = deli; //Point R4 at deli[0][0]
	P4.H = deli; //as usual, needs two ops
	P3.L = waco; //Point P3 at waco[0][0]
	P3.H = waco;
	
	R0 = 0; //Just for initializing stuff

	
	LSETUP(SETUP_LOOP, END_SETUP_LOOP) LC0=P2;
	
SETUP_LOOP:

//load up the frame with rstructs[], recap
/***
{
		[totalsamps						4 bytes, how many samples captured since start, level 0 only]
		insampcnt						4 bytes, counts how many samps fed in from level above waiting to be processed, should never exceed 2
		drr0 - address of deli[r]; 		4 bytes, origin of circular buffer we are convoluting 
		i0r - current I0[r];			4 bytes, where to write next two samples/sliding pointer in deli[r]
		wcr0 - address of waco[r];		4 bytes, origin of circular buffer for high-pass output, if we want to process on the fly
		i3r - current I3[r];			4 bytes, sliding pointer in waco[r], where we save high-pass output samples/wavelet coeffs
}
***/
	[SP--] = R0; //Zero the in-sample count for this level
	[SP--] = P4; //address of the deli[r][0]
	[SP--] = P4; //Initial position of I0 for this deli - at deli[r][0]
	P4 = P4 + P5; //advance P4 to point at next deli[]
	[SP--] = P3;	//address of waco[r][0]
	[SP--] = P3; 	//Initial position of I3 for waco circular buffer 
	P3 = P3 + P5;	//advance to next waco[]

END_SETUP_LOOP:	NOP;

//SP now points at the bottom of the frame containing the rstructs for each level, starting with level 0
	

//Get the ball rolling
	call Init_EBIU;
	call Init_Flash;
	call Init_1836;
	call Init_Sport0;
	call Init_DMA;
	call Init_Interrupts;
	call Enable_DMA_Sport0;
//Hereafter the ISR will be called every time the audio ADC DMA's us a sample (frame)

		
	R2 = [I2++]; //Get hh[0], hh[1] into R2.L, R2.H for the first time; thereafter they'll always be here after the completion of any hp convolution
	R3 = [I3++]; //Get hl[0], hl[1] into R3.L, R3.H  "   "    "     "   ........ " ......                                          "  lp    "
		
FOREVER:
/***********
No longer using 2-bit counters in a register and EXTRACT (but remember it for bit-wise operations). Addresses/counters for each level are in "structs"
in the frame
	//Get address of dd[1][] into P0
	R4.L = deli;
	R4.H = deli;
	R3 = TAPLENGTH<<1;
	R4 = R4 + R3;
	P0 = R4; //TAPLENGTH*2; //2-byte samples
	B0 = P0;
	I0 = P0;

****/
/*** Regarding filter indexing: use I2 for low-pass, I3 for high-pass, then once started both will run full circle in each convolution,
so after initial setup they don't need adjusting at any level. Caveat: the ISR uses them too, and may interrupt any lower-level cycle at any point.
Possible solution: save I2 and I3 to the (supervisor) stack, then reinitialize them to point at hl[0], hh[0] in every ISR invocation. No more costly
(I think) than pushing and popping them on the supervisor stack.
***/

	R4 = -2;
	P5 = SP;
	//Note that one cannot (apparently) subtract from address registers, except SP
	
//Check two-bit counter for loop 0. This is written to by the ISR, so must protect before accessing
	CLI R6;	//Disable interrupts
	R5 = [P5]; //Get the 2-bit counter for loop 0
	CC = R5 < 2; //Have (at least) 2 samples to go? CC will be false if we have
	IF CC JUMP NOTGOT2; //If not got 2 samples to work with don't clear R5. The CC will ensure no level-0 processing is done
	R5 = R5 + R4; // R4; // Get here if we have 2 samples to work on. 0 count, it will be written back
NOTGOT2:
	[P5++] = R5; //Write back (changed or not) and advance pointer in structs
//	STI R6; //Reneable interrupts - No, not yet, if level 0 is to be done it must complete before a new ISR can add samples and move I0

START_RLOOP:	
	LSETUP(RLOOP, END_RLOOP) LC0 = P2;//that's R-1, no. of resolution levels except the first, which is done [in the ISR?] now in userland before RLOOP

RLOOP:
//First check 2-bit counter to see if we have 2 samples to process. Note that LC0 will count down, can use it to refer to ...
//If this is loop 1, the counter can be incremented by the capture ISR, so we must avoid conflict by disabling interrupts  
	

//CC is flagging the outcome of testing in the previous level this level's two-bit counter, to see if we have two new samples yet
	IF !CC JUMP DO_TWO; //If got two samples to work with, do convolutions for this level
	P5 += HOP_OVER;	//Else, advance P5 to point at rstruct for the next resolution level, move to next (in dwt cascade lower resolution) level
	JUMP END_RLOOP;
	
//Do high-pass filter convolution first, yielding detail wavelet coeffs. The only reason is that after low-pass the (decimated) output sample can be
//immediately inserted into the delay-line of the next lower level
	
DO_TWO:	
//First reestablish the delay line for this level as a circular buffer. P5 is now pointing at the saved origin drr0 for level r
	P0 = [P5++]; //Origin of the deli[r]
	B0 = P0;	 //make this a circular buffer, L0 is always TAPLENGTH
	P1 = [P5++]; //Get the stored I0 index which is currently pointing at the byte after the two most recently appended samples, call then x[n], x[n+1]
	I0 = P1;
	I0 -= M0;	//Now pointing at x[n+1]
	R7 = I0;	//Cache this I0, need for LP filter as well
	R0.H = W[I0--]; //Get x[n+1] into R0.H, I0 now pointing at x[n]
	R0.L = W[I0--];	//x[n] into R0.L, I0 now pointing at x[n-1]
	
//x[n], x[n+1] have to be in R0.L, R0.H by the time we get here, 
//h[0] must be in R2.L, h[1] in R2.H, should be there since the last convolution
	
//m = 0, overture to INLOOP			
	A0 = R2.L*R0.L, 	//h[0].x[n]
	A1 = R2.L*R0.H ||	//h[0].x[n+1]
//legacy, from the streaming FIR filter	[I0++] = R0;//write the two new samples into d[], at the address indexed by I0, advance I0
	R0.H = W[I0--];		//R0.H = x[n-1], R0.L = x[n], decrement I1 by one 2-byte step
	
	LSETUP(INLOOP,END_INLOOP) LC1=P3>>1; //deal with 2 samples in every loop, 1 cycle per sample, as in the original fir example	

INLOOP:	
	A0 += R2.H*R0.H,	//add h[1].x[n-1] to A0	N.B. {Two MAC ops in one instruction must both use the same two registers, }
	A1 += R2.H*R0.L ||	//add h[1].x[n] to A1. 		 { e.g if R0 was R1 in the second instruction it would be illegal}
	R2 = [I2++] ||		//get h[2], h[3] into R2, so R2.L = h[2], R2.H = h[3]
	R0.L = W[I0--];//I1--];		//R0.L = x[n-2], R0.H = x[n-1];								
	
END_INLOOP:
	A0 += R2.L*R0.L,	//add h[2].x[n-2] to A0
	A1 += R2.L*R0.H ||	//add h[2].x[n-1] to A1
	R0.H = W[I0--]; //I1--];		//R0.L= x[n-2], R0.H=x[n-3]
	
//Rework needed, should R6 be used here???
//	R6.L = (A0 += R2.H*R0.H),	//add h[7].x[n-7] to A0, assuming all done copy low word to R6.L 
//	R6.H = (A1 += R2.H*R0.L) ||	//add h[7].x[n-6] to A1,      "     "   "  copy low word to R6.H
	A0 += R2.H*R0.H,
	A1 += R2.H*R0.L ||
	R2 = [I2++]; 				//R2.L = h[0], R2.H = h[1], for next r-loop
	
	I0 = R7; //I0 pointing at x[n+1] again, do it now because of "Dagreg read after write ... requires 4 extra cycles"
	A0 += A1;		//add and
	A0 = A0>>1 ||		//average (my decimation)
//Write low word of A0 to the detail (waco) circular buffer:
	P0 = [P5++]; //Origin of the waco buffer
	B3 = P0;	//Reinstate this waco buffer as circular buffer 3
	P1 = [P5]; //Getthe stored I3 index pointing at next location to write in waco buff
	I3 = P1;
	R0 = A0.W;  //This copies low 16 bits into R0.L, but note that R0.L = A0.W gives a compiler error, don't know why yet
	W [I3++] = R0.L;
	R1 = I3;	
	[P5++] = R1; //Put new I3 back in storage, advance the frame read pointer
	



/************************
obsolete code	
	A0 += A1;		//Average of pair
	A0 = A0 >>> 1;  //Divide sum by 2
	R6 = A0;		//Stash the approximation coeff in R6
//	R6 <<= 16;// PACK(R6.L, R6.L); //to get low word into high word

	CC = BITTST(R4, DONE_LP);
	IF CC JUMP DONE_BOTH; //If not done HP go do it
	BITSET(R4, DONE_LP); 		//else set flag to show we did the LP
	CLI R5; //Never set flags without disabling interrupts
	W [P0] = R4.L; //Write back bit flags
	STI R5;
end of obsolete code*************/

/********************************************
* Switch to the LP filter					*
*********************************************/


//Already have I0 pointing at x[n+1] again, get x[n], x[n+1] into R0, slide I0 back to point at x[n-1]
//...
	R0.H = W[I0--]; //Get x[n+1] into R0.H, I0 now pointing at x[n]
	R0.L = W[I0--];	//x[n] into R0.L, I0 now pointing at x[n-1]
//m = 0, overture to INLOOP_LP			
	A0 = R3.L*R0.L, 	//hl[0].x[n]
	A1 = R3.L*R0.H ||	//hl[0].x[n+1]
	R0.H = W[I0--];		//R0.H = x[n-1], R0.L = x[n], decrement I1 by one 2-byte step

	LSETUP(INLOOP_LP,END_INLOOP_LP) LC1=P3>>1; //deal with 2 samples in every loop, 1 cycle per sample, as in the original fir example	

INLOOP_LP:	
	A0 += R3.H*R0.H,	//add h[1].x[n-1] to A0	N.B. {Two MAC ops in one instruction must both use the same two registers, }
	A1 += R3.H*R0.L ||	//add h[1].x[n] to A1. 		 { e.g if R0 was R1 in the second instruction it would be illegal}
	R2 = [I3++] ||		//get h[2], h[3] into R2, so R2.L = h[2], R2.H = h[3]
	R0.L = W[I0--];//I1--];		//R0.L = x[n-2], R0.H = x[n-1];								
	
END_INLOOP_LP:
	A0 += R3.L*R0.L,	//add h[2].x[n-2] to A0
	A1 += R3.L*R0.H ||	//add h[2].x[n-1] to A1
	R0.H = W[I0--];		//R0.L= x[n-2], R0.H=x[n-3], etc.
	
		
	STI R6; //Reenable interrupts. This has no effect after level 0, interrupts only needed to be disabled during level 0, 
			//because the ISR can change I0 at that level
	
//Rework needed, should R6 be used here???
	R6.L = (A0 += R3.H*R0.H),	//add h[7].x[n-7] to A0, assuming all done copy low word to R6.L 
	R6.H = (A1 += R3.H*R0.L) ||	//add h[7].x[n-6] to A1,      "     "   "  copy low word to R6.H
	R2 = [I2++]; // ||				//R2.L = h[0], R2.H = h[1], for next r-loop

END_RLOOP:  //Check whether next level has got 2 more samples to process, will be indicated by CC
			//Hopefully nothing else happens here except continue or exit loop depending on LC0
	R5 = [P5++]; //Get two-bit counter
	CC = R5 < 2; //Note that > cannot be used here - funny quirks of this assembler!
	

//Get origin of next lower deli, and the save I0 where the decimated (average in my version) x[<something>] needs to be written; increment its sample cnt
	

DONE_BOTH:						
	CLI R5;
	R4.L= W[P0]; 
	CC = BITTST(R4, BYE); 
	IF !CC JUMP FOREVER (BP);
						