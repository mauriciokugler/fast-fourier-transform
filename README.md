# ANSI C++ Fast Fourier Transform

Yet another C++ implementation of the [Fast Fourier Transform](https://en.wikipedia.org/wiki/Fast_Fourier_transform) (FFT) algorithm. This ANSI C++ class provides simple and efficient methods for 1D, 2D & 3D direct and inverse FFT calculation. 

## Usage

The 1D-FFT (*n1* points), 2D-FFT (*n1*&#215;*n2* points) and 3D-FFT (*n1*&#215;*n2*&#215;*n3* points) objects can be instantiated with the following constructors:

```C++
FastFourierTransform(unsigned int n1)
FastFourierTransform(unsigned int n1, unsigned int n2)
FastFourierTransform(unsigned int n1, unsigned int n2, unsigned int n3)
```

The number of points in each dimension, *n1*, *n2* and *n3*, should be (possibly different) powers of 2. If any other value is given, it will be set to the next nearest power of 2. If any of these parameters is zero, the object won't be initialized. 

The direct Fourier Transform is calculated using the following methods:

```C++
complex<float> inline *fft1(float *x)
complex<float> inline **fft2(float **x)
complex<float> inline ***fft3(float ***x) 
```

The input array *x* should have the same dimensions of the parameters used in the constructor. For efficiency reasons, the values in the return array are the complex-conjugate of the actual values.

The inverse Fourier Transform is calculated using the following methods:

```C++
float inline *ifft1(complex<float> *x)  
float inline **ifft2(complex<float> **x) 
float inline ***ifft3(complex<float> ***x)  
```

The destructor method delete all allocated arrays, including the ones which pointers have been returned by the `fft` and `ifft` methods.  


## Example

The `main.cpp` file contains a simple example of the 2D-FFT:

```C++
#include "FastFourierTransform.h"
#include <omp.h>

int main()
{
  const unsigned int n1 = 8;
  const unsigned int n2 = 4;

  float **x = new float*[n1];
  for(unsigned int i=0;i<n1;i++) {
    x[i] = new float[n2];
    for(unsigned int j=0;j<n2;j++) {
      x[i][j] = ((float)rand()/(RAND_MAX))*2-1;
    }
  }

  FastFourierTransform *FFT = new FastFourierTransform(n1, n2);

  complex<float> **y = FFT->fft2(x);
  float **z = FFT->ifft2(y);

  delete FFT;

  return 0;
}
```

The [OpenMP](https://en.wikipedia.org/wiki/OpenMP) parallelization is optional, but has critical impact on the code's performance. 

## Citing

```TeX
@MISC{Kugler2019,
  author = "Mauricio Kugler",
  title = "ANSI C++ Fast Fourier Transform",
  year = "2019",
  url = "https://github.com/mauriciokugler/fast-fourier-transform",
  note = "Version 1.0.0"
}
```

## License

This project is licensed under the [MIT License](LICENSE).