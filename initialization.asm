#include "Talkthrough.h"
						
/*****************************************************************************
 Function:	Init_EBIU													
 
 Description:	This function initialises and enables asynchronous memory 	
		banks in External Bus Interface Unit so that Flash A can be
		accessed.												
******************************************************************************/
.section L1_data_a; //Introduced here by Paul for dwt_cascade. In ASM_Talkthrough example and ASM_Talkthrough_TDM_nach_Paul it was
					//included in main.asm, seemed more logical to have it here

// left input data from ad1836
.var Channel0LeftIn, Channel1LeftIn;
// right input data from ad1836
.var Channel0RightIn, Channel1RightIn;
// left ouput data for ad1836
.var Channel0LeftOut, Channel1LeftOut;
// right ouput data for ad1836
.var Channel0RightOut, Channel1RightOut;
// array for registers to configure the ad1836
// names are defined in "Talkthrough.h"
.align 2;
.byte2	Codec1836TxRegs[CODEC_1836_REGS_LENGTH] =
{									
					DAC_CONTROL_1	| 0x000,
					DAC_CONTROL_2	| 0x000,
					DAC_VOLUME_0	| 0x3ff,
					DAC_VOLUME_1	| 0x3ff,
					DAC_VOLUME_2	| 0x3ff,
					DAC_VOLUME_3	| 0x3ff,
					DAC_VOLUME_4	| 0x3ff,
					DAC_VOLUME_5	| 0x3ff,
					ADC_CONTROL_1	| 0x000,
					ADC_CONTROL_2	| 0x180,
					ADC_CONTROL_3	| 0x000
					
};


// SPORT0 DMA transmit buffer
.align 4;
.byte4	tx_buf[8];
// SPORT0 DMA receive buffer
.align 4;
.byte4	rx_buf[8];


.section L1_code;
Init_EBIU:

	P0.L = LO(EBIU_AMBCTL0);
	P0.H = HI(EBIU_AMBCTL0);
	R0.L = 0x7bb0;
	R0.H = 0x7bb0;
	[ P0 ] = R0;
	
	P0.L = LO(EBIU_AMBCTL1);
	P0.H = HI(EBIU_AMBCTL1);
	R0.L = 0x7bb0;
	R0.H = 0x7bb0;
	[ P0 ] = R0;
	
	P0.L = LO(EBIU_AMGCTL);
	P0.H = HI(EBIU_AMGCTL);
	R0 = 0x000f;
	W[ P0 ] = R0.L;

Init_EBIU.END:
	RTS;

/*****************************************************************************
 Function:	Init_Flash													
 
 Description:	This function initialises pin direction of Port A in Flash A
		to output.  The AD1836_RESET on the ADSP-BF533 EZ-KIT board 
		is connected to Port A.										
******************************************************************************/
.section L1_code;
Init_Flash:

	P0.L = LO(FlashA_PortA_Dir);
	P0.H = HI(FlashA_PortA_Dir);
	R0 = 0x1(z);
	W[ P0 ] = R0;
	
Init_Flash.END:
	RTS;

/*****************************************************************************
 Function:	Init1836													
 
 Description:	This function sets up the SPI port to configure the AD1836. 
		The content of the array Codec1836TxRegs is sent to the 
		codec.														
******************************************************************************/
.section L1_code;
Init_1836:

	// write to Port A to reset AD1836
	P0.L = LO(FlashA_PortA_Data);
	P0.H = HI(FlashA_PortA_Data);
	R0 = 0x00(z);
	B[ P0 ] = R0;

	// write to Port A to enable AD1836
	R0 = 0x01(z);
	B[ P0 ] = R0;

	// wait to recover from reset
	P5.H = 0;
	P5.L = 0xfff0;
	loop reset_recover LC0 = P5;
	LOOP_BEGIN reset_recover;
		nop;
		nop;
		nop;
		nop;
		nop;
		nop;
	LOOP_END reset_recover;

	// Enable PF4
	P0.L = LO(SPI_FLG);
	P0.H = HI(SPI_FLG);
	R0.L = W[ P0 ];
	BITSET(R0, FLS4_P);
	W[ P0 ] = R0;
	
	// Set baud rate SCK = HCLK/(2*SPIBAUD)
	P0.L = LO(SPI_BAUD);
	R0 = 16;
	W[ P0 ] = R0.L;
	
	// configure spi port
	// SPI DMA write, 16-bit data, MSB first, SPI Master
	P0.L = LO(SPI_CTL);
	R0 = TIMOD_DMA_TX | SIZE | MSTR;
	W[ P0 ] = R0.L;
	
	// Set up DMA5 to transmit
	// Map DMA5 to SPI
	P1.L = LO(DMA5_PERIPHERAL_MAP);
	P1.H = HI(DMA5_PERIPHERAL_MAP);
	R1.L = 0x5000;
	W[ P1 ] = R1;
	
	// Configure DMA5
	// 16-bit transfers
	P1.L = LO(DMA5_CONFIG);
	R1.L = WDSIZE_16;			
	W[ P1 ] = R1.L; 
  	
	P2.H = Codec1836TxRegs;
	P2.L = Codec1836TxRegs;
	R1 = P2;
	
	// Start address of data buffer
	P1.L = LO(DMA5_START_ADDR);
	[ P1 ] = R1; 

	// DMA inner loop count
	P1.L = LO(DMA5_X_COUNT);
	R1.L = CODEC_1836_REGS_LENGTH;
	W[ P1 ] = R1.L; 

	// Inner loop address increment
	P1.L = LO(DMA5_X_MODIFY);
	R1.L = 2;
	W[ P1 ] = R1.L; 

	// Enable DMA
	P1.L = LO(DMA5_CONFIG);
	R1.L = W[ P1 ];
	BITSET(R1,DMAEN_P);
	W[ P1 ] = R1.L; 

	// ENABLE SPI
	P1.H = HI(SPI_CTL);
	P1.L = LO(SPI_CTL);
	R1.L = W[ P1 ];
	BITSET(R1,14);
	W[ P1 ] = R1; 
	
	// wait until DMA transfers for spi are finished
	P5.L = 0xf000;
	loop dma_finish LC0 = P5;
	loop_begin dma_finish;
		nop;
		nop;
		nop;
		nop;
		nop;
		nop;
	loop_end dma_finish;

	// DISABLE SPI
	P1.H = HI(SPI_CTL);
	P1.L = LO(SPI_CTL);
	R1.L = 0;
	W[ P1 ] = R1; 
	
Init_1836.END:
	RTS;


/*****************************************************************************
 Function:	Init_Sport0													
 
 Description:	Configure Sport0 for TDM mode, to transmit/receive data 	
		to/from the AD1836. Configure Sport for external clocks and 
		framesyncs.												
******************************************************************************/
.section L1_code;
Init_Sport0:	

	P0.H = HI(SPORT0_RCR1);
		
	// Sport0 receive configuration
	// External CLK, External Frame sync, MSB first
	// 32-bit data
	P0.L = LO(SPORT0_RCR1);
	R0 = RFSR;
	W[ P0 ] = R0.L;
	
	P0.L = LO(SPORT0_RCR2);
	R0 = SLEN_32;
	W[ P0 ] = R0.L;

	// Sport0 transmit configuration
	// External CLK, External Frame sync, MSB first
	// 24-bit data
	P0.L = LO(SPORT0_TCR1);
	R0 = TFSR;
	W[ P0 ] = R0.L;
	
	P0.L = LO(SPORT0_TCR2);
	R0 = SLEN_32;
	W[ P0 ] = R0.L;
	
	// Enable MCM 8 transmit & receive channels
	P0.L = LO(SPORT0_MTCS0);
	R0 = 0x000000FF;
	[ P0 ] = R0;
	
	P0.L = LO(SPORT0_MRCS0);
	R0 = 0x000000FF;
	[ P0 ] = R0;
	
	// Set MCM configuration register and enable MCM mode
	P0.L = LO(SPORT0_MCMC1);
	R0 = 0x0000;
	W[ P0 ] = R0;
	
	P0.L = LO(SPORT0_MCMC2);
	R0 = 0x101c;
	W[ P0 ] = R0;
	
Init_Sport0.END:
	RTS;
	
	
/*****************************************************************************
 Function:	Init_DMA												
 
 Description:	Initialise DMA1 in autobuffer mode to receive and DMA2 in	
				autobuffer mode to transmit									
******************************************************************************/
.section L1_code;
Init_DMA:
	
	// Set up DMA1 to receive
	// Map DMA1 to Sport0 RX
	P1.L = LO(DMA1_PERIPHERAL_MAP);
	P1.H = HI(DMA1_PERIPHERAL_MAP);
	R1.L = 0x1000;
	W[ P1 ] = R1;
	
	// Configure DMA1
	// 32-bit transfers, Interrupt on completion, Autobuffer mode
	P1.L = LO(DMA1_CONFIG);
	R1 = WNR | WDSIZE_32 | DI_EN | FLOW_1;	
	W[ P1 ] = R1.L; 
  	
	P2.H = rx_buf;
	P2.L = rx_buf;
	R1 = P2;
	
	// Start address of data buffer
	P1.L = LO(DMA1_START_ADDR);
	[ P1 ] = R1; //Puts the DMA dest addr in to [P1]

	// DMA inner loop count
	P1.L = LO(DMA1_X_COUNT);
	R1.L = 8;	
	W[ P1 ] = R1.L; 

	// Inner loop address increment
	P1.L = LO(DMA1_X_MODIFY);
	R1.L = 4;
	W[ P1 ] = R1.L; 
	
	// Set up DMA2 to transmit
	// Map DMA2 to Sport0 TX
	P1.L = LO(DMA2_PERIPHERAL_MAP);
	P1.H = HI(DMA2_PERIPHERAL_MAP);
	R1.L = 0x2000;
	W[ P1 ] = R1;
	
	// Configure DMA2
	// 32-bit transfers, Autobuffer mode
	P1.L = LO(DMA2_CONFIG);
	R1.L = WDSIZE_32 | FLOW_1;			
	W[ P1 ] = R1.L; 
  	
	P2.H = tx_buf;
	P2.L = tx_buf;
	R1 = P2;
	
	// Start address of data buffer
	P1.L = LO(DMA2_START_ADDR);
	[ P1 ] = R1; 

	// DMA inner loop count
	P1.L = LO(DMA2_X_COUNT);
	R1.L = 8;
	W[ P1 ] = R1.L; 

	// Inner loop address increment
	P1.L = LO(DMA2_X_MODIFY);
	R1.L = 4;
	W[ P1 ] = R1.L; 

Init_DMA.END:
	RTS;

	
/*****************************************************************************
 Function:	Init_Interrupts											
 
 Description:	Initialise Interrupt for Sport0 RX							
******************************************************************************/
.section L1_code;
Init_Interrupts:

	[--SP] = R7; //apparently Push Multiple (q.v.) doesn't work with P regs on the BF533
	[--SP] = R1;
	[--SP] = R0;
	[--SP] = P1;
	[--SP] = P0;

	// Set Sport0 RX (DMA1) interrupt priority to 2 = IVG9 
	P0.L = LO(SIC_IAR1);
	P0.H = HI(SIC_IAR1);	
	R1.L = 0xff2f;
	R1.H = 0xffff;
	[ P0 ] = R1;
	
	// Unmask peripheral SPORT0 RX interrupt
	P0.L = LO(SIC_IMASK);
	P0.H = HI(SIC_IMASK);
	R1 = [ P0 ];
	BITSET(R1, 9);
	[ P0 ] = R1;

	// Remap the vector table pointer from the default __I9HANDLER 
	// to the new _SPORT0_RX_ISR interrupt service routine
	P0.L = LO(EVT9);
	P0.H = HI(EVT9);
	R0.l = _SPORT0_RX_ISR;
	R0.h = _SPORT0_RX_ISR;
	[ P0 ] = R0;
	
	// Enable interrupts IVG9
	P0.L = LO(IMASK);
	P0.H = HI(IMASK);
	R7 = [ P0 ];
	R1.H = 0;
	R1.L = 0x0200;
	R7 = R7 | R1;
	[ P0 ] = R7;
	
	P0 = [SP++];
	P1 = [SP++];
	R0 = [SP++];
	R1 = [SP++];
	R7 = [SP++];
	
	
Init_Interrupts.END:
	RTS;
	
/*****************************************************************************
 Function:	Enable_DMA_Sport										
 
 Description:	Enable DMA1, DMA2, Sport0 TX and Sport0 RX					
******************************************************************************/
.section L1_code;
Enable_DMA_Sport0:

	[--SP] = P1; 
	[--SP] = R1;

	// Enable DMA2
	P1.L = LO(DMA2_CONFIG);
	R1.L = W[ P1 ];
	BITSET(R1,DMAEN_P);
	W[ P1 ] = R1.L; 

	// Enable DMA1
	P1.L = LO(DMA1_CONFIG);
	R1.L = W[ P1 ];
	BITSET(R1,DMAEN_P);
	W[ P1 ] = R1.L; 
	
	// ENABLE SPORT0 TX
	P1.H = HI(SPORT0_TCR1);
	P1.L = LO(SPORT0_TCR1);
	R1.L = W[ P1 ];
	BITSET(R1,0);
	W[ P1 ] = R1; 

	// ENABLE SPORT0 RX
	P1.L = LO(SPORT0_RCR1);
	R1.L = W[ P1 ];
	BITSET(R1,0);
	W[ P1 ] = R1; 
	
	R1 = [SP++];
	P1 = [SP++];
	
Enable_DMA_Sport0.END:
	RTS;
	
	
