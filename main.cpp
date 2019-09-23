//Copyright (c) 2019 Mauricio Kugler, Nagoya Institute of Technology

#include "FastFourierTransform.h"
#include <stdlib.h>
#include <omp.h>

int main()
{
	const unsigned int n1 = 8;
	const unsigned int n2 = 4;

	float **x = new float*[n1];
	for(unsigned int i=0;i<n1;i++) {
		x[i] = new float[n2];
		for(unsigned int j=0;j<n2;j++) {
			x[i][j] = ((float)rand()/RAND_MAX)*2-1;
		}
	}

	FastFourierTransform *FFT = new FastFourierTransform(n1, n2);

	complex<float> **y = FFT->fft2(x);
	float **z = FFT->ifft2(y);

	delete FFT;

	return 0;
}