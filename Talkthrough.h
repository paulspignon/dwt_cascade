//--------------------------------------------------------------------------//
// Header files	
// Imported to cascade_dwt with (as yet) mininal changes from the Talkthrough example															//
//--------------------------------------------------------------------------//
#include <defBF533.h>

//--------------------------------------------------------------------------//
// Symbolic constants														//
//--------------------------------------------------------------------------//
// addresses for Port B in Flash A
#define FlashA_PortA_Dir	0x20270006
#define FlashA_PortA_Data	0x20270004

// names for codec registers, used for Codec1836TxRegs[]
#define DAC_CONTROL_1		0x0000
#define DAC_CONTROL_2		0x1000
#define DAC_VOLUME_0		0x2000
#define DAC_VOLUME_1		0x3000
#define DAC_VOLUME_2		0x4000
#define DAC_VOLUME_3		0x5000
#define DAC_VOLUME_4		0x6000
#define DAC_VOLUME_5		0x7000
#define ADC_0_PEAK_LEVEL	0x8000
#define ADC_1_PEAK_LEVEL	0x9000
#define ADC_2_PEAK_LEVEL	0xA000
#define ADC_3_PEAK_LEVEL	0xB000
#define ADC_CONTROL_1		0xC000
#define ADC_CONTROL_2		0xD000
#define ADC_CONTROL_3		0xE000

// names for slots in ad1836 audio frame
#define INTERNAL_ADC_L0			0
#define INTERNAL_ADC_L1			4
#define INTERNAL_ADC_R0			16
#define INTERNAL_ADC_R1			20
#define INTERNAL_DAC_L0			0
#define INTERNAL_DAC_L1			4
#define INTERNAL_DAC_L2			8
#define INTERNAL_DAC_R0			16
#define INTERNAL_DAC_R1			20
#define INTERNAL_DAC_R2			24


// size of array Codec1836TxRegs
#define CODEC_1836_REGS_LENGTH	11

// SPI transfer mode
#define TIMOD_DMA_TX 0x0003

// SPORT0 word length
#define SLEN_32	0x001f

// DMA flow mode
#define FLOW_1	0x1000


//--------------------------------------------------------------------------//
// Global variables															//
//--------------------------------------------------------------------------//
.extern Channel0LeftIn;
.extern Channel0RightIn;
.extern Channel0LeftOut;
.extern Channel0RightOut;
.extern Channel1LeftIn;
.extern Channel1RightIn;
.extern Channel1LeftOut;
.extern Channel1RightOut;
.extern Codec1836TxRegs;
.extern rx_buf;
.extern tx_buf;


//--------------------------------------------------------------------------//
// Prototypes																//
//--------------------------------------------------------------------------//
// in file Initialisation.asm
.extern Init_EBIU;
.extern Init_Flash;
.extern Init_1836;
.extern Init_Sport0;
.extern Init_DMA;
.extern Init_Interrupts;
.extern Enable_DMA_Sport0;

// in file interrupts.asm
.extern _SPORT0_RX_ISR;

// in file Process_data.asm
//.extern Process_Data;


