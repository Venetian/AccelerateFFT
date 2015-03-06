//
//  accFFT.cpp
//  AccelerateFFTtool
//
//  Created by Andrew Robertson on 17/07/2012.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "accFFT.h"


accFFT :: accFFT(int fft_size,int type)
{
    fft_type = type;//see below, 0 is float, 1 is double
    
    fftSize = fft_size;           
    fftSizeOver2 = fftSize/2;
    log2n = log2f(fftSize);  
    log2nOver2 = log2n/2;
    
    if (fft_type == 0)
    {
        split.realp = (float *) malloc(fftSize * sizeof(float));
        split.imagp = (float *) malloc(fftSize * sizeof(float));
        
        // allocate the fft object once
        fftSetup = vDSP_create_fftsetup(log2n, FFT_RADIX2);
        if (fftSetup == 0) {//was NULL
            printf("FFT Setup failed\n");
        }
    }
    else if (fft_type == 1) 
    {
        d_split.realp = (double *) malloc(fftSize * sizeof(double));
        d_split.imagp = (double *) malloc(fftSize * sizeof(double));
        
        // allocate the fft object once
        fftSetupD = vDSP_create_fftsetupD(log2n, FFT_RADIX2);
        if (fftSetupD == 0){//NULL) {
            printf("FFT Setup failed\n");
        }
    }
        
    
    
}

accFFT :: ~accFFT()
{
    if (fft_type == 0)
    {
        free(split.realp);
        free(split.imagp);
        vDSP_destroy_fftsetup(fftSetup);
    }
    else if (fft_type == 1)
    {
        free(d_split.realp);
        free(d_split.imagp);
        vDSP_destroy_fftsetupD(fftSetupD);
    }
    
    
}

void accFFT :: forward_FFT_f(float *buffer,float *real,float *imag)
{        
    //convert to split complex format with evens in real and odds in imag
    vDSP_ctoz((COMPLEX *) buffer, 2, &split, 1, fftSizeOver2);
    //this is basically a format required for the FFT in apple's DSP 
	//so we cast buffer (real-valued) as a complex, then create the split version
	
	
    //calc fft
    vDSP_fft_zrip(fftSetup, &split, 1, log2n, FFT_FORWARD);
    
    // set Nyquist component to imaginary of 0 component
    split.realp[fftSizeOver2] = split.imagp[0];
    split.imagp[fftSizeOver2] = 0.0;
    
    // set 0 component to zero
    split.imagp[0] = 0.0;
    
    // multiply by 0.5 to get correct output (to do with Apple's FFT implementation)
    for (i = 0; i <= fftSizeOver2; i++)
    {
        split.realp[i] *= 0.5;
        split.imagp[i] *= 0.5;
    }
    
    // set values above N/2+1 which are complex conjugate mirror image of those below
    for (i = fftSizeOver2 - 1;i > 0;--i)
    {
        split.realp[2*fftSizeOver2 - i] = split.realp[i];
        split.imagp[2*fftSizeOver2 - i] = -1*split.imagp[i];
        
        //cout << split_data.realp[2*fftSizeOver2 - i] << "   " << split_data.imagp[2*fftSizeOver2 - i] << "i" << endl;
    }
    
    for (i = 0;i < fftSize;i++)
    {
        real[i] = split.realp[i];
        imag[i] = split.imagp[i];
    }
}



void accFFT :: forward_FFT_d(double *buffer,fft_complex *out)
{        
    //convert to split complex format with evens in real and odds in imag
    vDSP_ctozD((DOUBLE_COMPLEX *) buffer, 2, &d_split, 1, fftSizeOver2);
    
    //calc fft
    vDSP_fft_zripD(fftSetupD, &d_split, 1, log2n, FFT_FORWARD);
    
    // set Nyquist component to imaginary of 0 component
    d_split.realp[fftSizeOver2] = d_split.imagp[0];
    d_split.imagp[fftSizeOver2] = 0.0;
    
    // set 0 component to zero
    d_split.imagp[0] = 0.0;
    
	// multiply by 0.5 to get correct output (to do with Apple's FFT implementation)
    double scale = 0.5;
	vDSP_vsmulD(d_split.realp, 1, &scale, d_split.realp, 1, fftSizeOver2);
	vDSP_vsmulD(d_split.imagp, 1, &scale, d_split.imagp, 1, fftSizeOver2);
	/*
	 //more obvious way
    for (i = 0; i <= fftSizeOver2; i++)
    {
        d_split.realp[i] *= 0.5;
        d_split.imagp[i] *= 0.5;
    }
    
	 */
    // set values above N/2+1 which are complex conjugate mirror image of those below
    for (i = fftSizeOver2 - 1;i > 0;--i)
    {
        d_split.realp[2*fftSizeOver2 - i] = d_split.realp[i];
        d_split.imagp[2*fftSizeOver2 - i] = -1*d_split.imagp[i];
    }
    
    for (i = 0;i < fftSize;i++)
    {
        out[i][0] = d_split.realp[i]; 
        out[i][1] = d_split.imagp[i];
    }
}

void accFFT::inverse_FFT_example_d(double *buffer){
	
	//convert to split complex format with evens in real and odds in imag
    vDSP_ctozD((DOUBLE_COMPLEX *) buffer, 2, &d_split, 1, fftSizeOver2);
	
	for (int i = 0; i < fftSize; i++){
		printf("double complex [%i] %f + %fi\n", i, d_split.realp[i], d_split.imagp[i]);
	}
    
    //calc fft
    vDSP_fft_zripD(fftSetupD, &d_split, 1, log2n, FFT_FORWARD);
	
	//CONVERT APPLE'S COMPACT REPRESENTATION INTO PROPER FFT
	//SEE Apple documentation on vDSP for explanation
	//effectively, the imag value at [0] is always 0, so stores the value at Nyquist (N/2) instead
	//firstly we need to scale by 0.5 when we have done complex forward fft
	//due to the way apple performs
	double scale = 0.5;
	vDSP_vsmulD(d_split.realp, 1, &scale, d_split.realp, 1, fftSizeOver2);
	vDSP_vsmulD(d_split.imagp, 1, &scale, d_split.imagp, 1, fftSizeOver2);
	
	for (int i = 0; i < fftSizeOver2; i++){
		printf("fft_split[%i] %f + %fi\n", i, d_split.realp[i], d_split.imagp[i]);
	}
	//get rid of nasty nyquist value
	d_split.realp[fftSizeOver2] = d_split.imagp[0];
	d_split.imagp[0] = 0;
	d_split.imagp[fftSizeOver2] = 0;
	//end CONVERT
	
	//FFT has been calculated above
	
	//option to add to spectrum here
	addToMagnitudeSpectrum();
	
	
	//here is the inverse FFT
	vDSP_fft_zripD(fftSetupD, &d_split, 1, log2n, FFT_INVERSE);
	
	//secondly, scale back by n when doing an n size fft for inverse
	scale = (float) 1.0 / (fftSize);
	vDSP_vsmulD(d_split.realp, 1, &scale, d_split.realp, 1, fftSizeOver2);
	vDSP_vsmulD(d_split.imagp, 1, &scale, d_split.imagp, 1, fftSizeOver2);
	
	for (int i = 0; i < fftSizeOver2; i++){
		printf("split[%i] %f + %fi\n", i, d_split.realp[i], d_split.imagp[i]);
	}
	 
	
}

void accFFT::inverse_FFT_f(float* buffer){
	vDSP_ctoz((COMPLEX *) buffer, 2, &split, 1, fftSizeOver2);
    //see above with double - To be completed
}

void accFFT::addToMagnitudeSpectrum(){
	float* magVector = new float[fftSizeOver2];
	for (int i = 0; i < fftSizeOver2; i++){
		printf("FFT[%i] real %f, imag %f\n", i, d_split.realp[i], d_split.imagp[i]);
		magVector[i] = sqrt((d_split.realp[i]*d_split.realp[i]) + (d_split.imagp[i]*d_split.imagp[i]));
		printf("mag spectrum [%i] %f\n", i, magVector[i]);
	}
	magSpectrum.push_back(magVector);
}