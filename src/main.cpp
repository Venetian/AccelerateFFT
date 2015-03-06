#include <iostream>



#include <stdio.h>

#include <stdlib.h>

#include <Carbon/Carbon.h>

/* Carbon #include <Types.h> */

/* #include <OSUtils.h> */

/* Carbon #include <Timer.h> */

/* Carbon #include <Gestalt.h> */

/* #include <Errors.h> */

/* #include <vfp.h> */

#include <Accelerate/Accelerate.h>

#include <sys/time.h>

#include <iostream>

/* #include "vDSP.h" */

#include "accFFT.h"
using namespace std;


#define MAX_LOOP_NUM       10000 /* Number of iterations used in

the timing loop */



/* used in looking for a g4 */

#define kHasAltiVecMask    ( 1 << gestaltPowerPCHasVectorInstructions  )



#define MAX(a,b)           ( (a>b)?a:b )

#define N                  2  /* This is a power of 2 defining  the length

* of the FFTs */

FFTSetup fftSetup;
COMPLEX_SPLIT FFT_A;



static struct timeval startTime, endTime;



#if (defined(__ppc__) || defined(__ppc64))



#define kHasAltiVecMask ( 1 << gestaltPowerPCHasVectorInstructions )



Boolean         HasAltivec();

void            TurnJavaModeOff(vector unsigned int *oldJavaMode);

void            RestoreJavaMode(vector unsigned int *oldJavaMode);



#endif



void            Start_Clock(void);

void            Stop_Clock(float *call_time);

void            Compare(float *original, float *computed, long length);

void            Dummy_fft_zip(FFTSetup setup, COMPLEX_SPLIT * C, long stride, unsigned long log2n, long flag);

void            Dummy_fft2d_zip(FFTSetup setup, COMPLEX_SPLIT * A, long rowStride, long columnStride, unsigned log2nr, unsigned log2nc, long flag);

void            Dummy_fft_zrip(FFTSetup setup, COMPLEX_SPLIT * C, long stride, unsigned long log2n, long flag);

void            Dummy_ctoz(COMPLEX * C, int32_t strideC, COMPLEX_SPLIT * Z, int32_t strideZ, int32_t size);

void            Dummy_ztoc(COMPLEX_SPLIT * A, int32_t strideZ, COMPLEX * C, int32_t strideC, int32_t size);

void            RealFFTUsageAndTiming();

void            ComplexFFTUsageAndTiming();

void            Complex2dFFTUsageAndTiming();





void doFFTReal(float samples[], float amp[], int numSamples)
{

	vDSP_Length log2n = log2f(numSamples);
	fftSetup = vDSP_create_fftsetup(log2n, FFT_RADIX2);
	int nOver2 = numSamples/2;
	FFT_A.realp = (float *) malloc(nOver2*sizeof(float));
	FFT_A.imagp = (float *) malloc(nOver2*sizeof(float));
	
	
	//from FFTaccelerate on github by Thomas Hobbes http://www.tomhoddes.com/?p=12
	int i;
	//vDSP_Length log2n = log2f(numSamples);
	
    //Convert float array of reals samples to COMPLEX_SPLIT array A
	vDSP_ctoz((COMPLEX*) samples, 2, &FFT_A, 1, numSamples/2);
	
    //Perform FFT using fftSetup and A
    //Results are returned in A
	vDSP_fft_zrip(fftSetup, &FFT_A, 1, log2n, FFT_FORWARD);
    
    //Convert COMPLEX_SPLIT A result to float array to be returned
    amp[0] = FFT_A.realp[0]/(numSamples*2);
	printf("Tom Hobbes[%i] : %f\n", 0, amp[0]);
	for(i=1;i<numSamples;i++){
		amp[i]=sqrt(FFT_A.realp[i]*FFT_A.realp[i]+FFT_A.imagp[i]*FFT_A.imagp[i])/numSamples;
		printf("Tom Hobbes[%i] : %f\n", i, amp[i]);
	}
}




int main()

{
	
    //RealFFTUsageAndTiming();
	
    /* ComplexFFTUsageAndTiming ( ); */
	
    /* Complex2dFFTUsageAndTiming ( ); */
	
	//////////////////////////////////////////
    accFFT fft(8,1);
    
    double my_buff[8] = {1,4,2,1,1,7,7,9};//this is alternating real, imaginary,...
//    double real[8];
  //  double complex[8];
    
	double complexResult[8][2];
    
	fft.inverse_FFT_example_d(my_buff);
	
    fft.forward_FFT_d(my_buff, complexResult);
    
    
    for (int i = 0;i < 8;i++)
    {
        cout << "real " << complexResult[i][0] << " imag  " << complexResult[i][1] << endl;
    }
    
	//////////////////////////////////////////    
	
	
	
	
    return 0;
	
}



/************************************************************************
 
 Real FFT.
 
 
 
 The real FFT, unlike the complex FFT, may possibly have to  use two
 
 transformation functions, one before the FFT call and one after.
 
 This is if the input array is not in the even-odd split configuration.
 
 
 
 A real array A = {A[0],...,A[n]} has to be transformed into  an
 
 even-odd array
 
 
 
 AEvenOdd = {A[0],A[2],...,A[n-1],A[1],A[3],...A[n]}
 
 
 
 via the vDSP_ctoz call.
 
 
 
 The result of the real FFT of AEvenOdd of dimension n is a  complex
 
 array of the dimension 2n, with a special format:
 
 
 
 {[DC,0],C[1],C[2],...,C[n/2],[NY,0],Cc[n/2],...,Cc[2],Cc[1]}
 
 
 
 where
 
 
 
 1. DC and NY are the dc and nyquist components (real valued),
 
 2. C is complex in a split representation,
 
 3. Cc is the complex conjugate of C in a split representation.
 
 
 
 For an n size real array A, the complex results require 2n
 
 spaces.  In order to fit the 2n size result into an n size  input and
 
 since the complex conjugates are duplicate information, the  real
 
 FFT produces its results as follows:
 
 
 
 {[DC,NY],C[1],C[2],...,C[n/2]}.
 
 
 
 The time for a length 1024 real FFT with the transformation  functions
 
 on a 500mhz machine is 14.9 microseconds.  If the data structure of the  calling
 
 program renders the transformation functions unecessary, the  same
 
 signal is processed in 11.3 microseconds.
 
 
 
 ************************************************************************/



void

RealFFTUsageAndTiming()

{
	
    COMPLEX_SPLIT   A;	
//	typedef DSPComplex                      COMPLEX;

	
    FFTSetup        setupReal;
	
    uint32_t        log2n;
	
    uint32_t        n, nOver2;
	
    int32_t         stride;
	
    uint32_t        i;
	
    float          *originalReal, *obtainedReal;
	
    float           scale;
	
	
	
    /* Set the size of FFT. */
	
    log2n = N;
	
    n = 1 << log2n;
	
	printf("log2n %i N %i\n", log2n, n);
	
    stride = 1;//process every element of row
	
    nOver2 = n / 2;
	
	
	
    printf("1D real FFT of length log2 ( %d ) = %d\n\n", n, log2n);
	
	
	
    /* Allocate memory for the input operands and check its availability,
	 
     * use the vector version to get 16-byte alignment. */
	
    A.realp = (float *) malloc(nOver2 * sizeof(float));
	
    A.imagp = (float *) malloc(nOver2 * sizeof(float));
	
    originalReal = (float *) malloc(n * sizeof(float));
	
    obtainedReal = (float *) malloc(n * sizeof(float));
	
	
	
    if (originalReal == NULL || A.realp == NULL || A.imagp == NULL) {
		
        printf("\nmalloc failed to allocate memory for  the real FFT"
			   
               "section of the sample.\n");
		
        exit(0);
		
    }
	
    /* Generate an input signal in the real domain. */
	
    for (i = 0; i < n; i++){
		
        originalReal[i] = 0;//(float) (i + 1);
		originalReal[0] = 1;
		//originalReal[2] = 0;
		printf("orig[%i] = %f\n", i, originalReal[i]);
	}

	//alt version via Tom Hobbes
	float fftOutput[n];
	//doFFTReal(originalReal, fftOutput, n);
	//end tom hobbes
	
    /* Look at the real signal as an interleaved complex vector  by
	 
     * casting it.  Then call the transformation function vDSP_ctoz to
	 
     * get a split complex vector, which for a real signal, divides into
	 
     * an even-odd configuration. */
	
    vDSP_ctoz((COMPLEX *) originalReal, 2, &A, 1, nOver2);
    //convert to split complex format with evens in real and odds in imag

	
/*	for (int i = 0;i < n;i++){
		printf("complex version[%i] = %f + %fi\n", i,  A.realp[i], A.imagp[i]);
	}
*/	
    /* Set up the required memory for the FFT routines and check  its
	 
     * availability. */
	
    setupReal = vDSP_create_fftsetup(log2n, FFT_RADIX2);
	
    if (setupReal == NULL) {
		
        printf("\nFFT_Setup failed to allocate enough memory  for"
			   
               "the real FFT.\n");
		
        exit(0);
		
    }
	for (int i = 0;i < n;i++){
		printf("input[%i] = %f + %fi\n", i,  A.realp[i], A.imagp[i]);
	}
    /* Carry out a Forward and Inverse FFT transform. */
	
	//MAIN FFT CALL
    vDSP_fft_zrip(setupReal, &A, stride, log2n, FFT_FORWARD);
	//AR NB use zripD for double
	
	//intermediate output
	for (int i = 0;i < n;i++){
		printf("output[%i] = %f + %fi\n", i,  A.realp[i], A.imagp[i]);
	}
	
	
	for (i = 0; i < nOver2; ++i)
    {
        A.realp[i] *= 0.5;
        A.imagp[i] *= 0.5;
    }
		
	printf("UNPACKED OUTPUT\n");
	printf("	(%d): %8g + %8g i\n", 0, A.realp[0], 0.0); // DC 0 or ?? A.imagp[0]

	for (i = 1; i < nOver2; ++i)
    {
        printf("(%d): %8g + %8g i\n", i, A.realp[i], A.imagp[i]);
    }
	
	
	//AR comment out inverse
  //  vDSP_fft_zrip(setupReal, &A, stride, log2n, FFT_INVERSE);
	
	
	
    /* Verify correctness of the results, but first scale it by  2n. */
	/*
    scale = (float) 1.0 / (2 * n);
	
    vDSP_vsmul(A.realp, 1, &scale, A.realp, 1, nOver2);
	
    vDSP_vsmul(A.imagp, 1, &scale, A.imagp, 1, nOver2);
	
	*/
	
    /* The output signal is now in a split real form.  Use the  function
	 
     * vDSP_ztoc to get a split real vector. */
	

    vDSP_ztoc(&A, 1, (COMPLEX *) obtainedReal, 2, nOver2);
	
	/*
	for (int i = 0;i < n;i++){
		printf("obtained real[%i] = %f ", i, obtainedReal[i], A.imagp[i]);
	}
	 */
	
    /* Check for accuracy by looking at the inverse transform  results. */
	
   // Compare(originalReal, obtainedReal, n);

	/* Free the allocated memory. */
	
    vDSP_destroy_fftsetup(setupReal);
	
    free(obtainedReal);
	
    free(originalReal);
	
    free(A.realp);
	
    free(A.imagp);
	
}




