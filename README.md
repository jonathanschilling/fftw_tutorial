# FFTW Tutorial

This is a basic C project (Makefile, but also for the Eclipse IDE) I use for exploring FFTW.

One- and two-dimensional DFTs of random data are computed using both FFTW and straight-forward naive algorithms
in order to illustrate explicitly what kind of symmetries and scaling properties FFTW implies in its inputs and outputs.

## One-Dimensional Examples
This tutorial starts by computing one-dimensional (1D) DFTs of random input data.

### 1D complex-to-complex
The first example is basically a self-contained version of the [corresponding example in the FFTW manual](http://fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html#Complex-One_002dDimensional-DFTs). It is available in the file [`test_1d_c2c.c`](https://github.com/jonathanschilling/fftw_tutorial/blob/master/src/test_1d_c2c.c).

```C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

void test_1d_c2c() {
    int n = 32;
    int dir = FFTW_BACKWARD;
    double real, imag, phi;

    fftw_complex *in = fftw_alloc_complex(n);
    fftw_complex *ref_out = fftw_alloc_complex(n);
    fftw_complex *fftw_out = fftw_alloc_complex(n);

    fftw_plan p = fftw_plan_dft_1d(n, in, fftw_out, dir, FFTW_ESTIMATE);

    // fill the input array with random data
    fill_random_1d_cplx(n, in);

    // compute the reference output
    for (int k = 0; k < n; ++k) {
        ref_out[k] = 0.0;
        for (int j = 0; j < n; ++j) {
            phi = dir * 2.0 * M_PI * j * k / ((double) n);

            real = creal(in[j]) * cos(phi) - cimag(in[j]) * sin(phi);
            imag = creal(in[j]) * sin(phi) + cimag(in[j]) * cos(phi);
            ref_out[k] += real + I * imag;
        }
    }

    // compute the DFT of in using FFTW
    fftw_execute(p);

    // compare reference output with FFTW output
    double eps = 1e-12;
    compare_1d_cplx(n, ref_out, fftw_out, eps);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(ref_out);
    fftw_free(fftw_out);
}

int main(int argc, char** argv) {
    test_1d_c2c();
    return 0;
}
```

## Allocation of arrays
Throughout this example collection, the proposed convenience wrapper functions provided by FFTW for allocating real- and complex-valued arrays are used:
```C
int n = 32;
int nOut = n/2+1;
double *in = fftw_alloc_real(n);
fftw_complex *out = fftw_alloc_complex(nOut);
```
where `N` is the real-space size of the DFT and `outN` is the number of Fourier coefficients resulting from a `r2c` DFT.
The corresponding "raw" allocation code would look like this:
```C
double *in = (double*) fftw_malloc(sizeof(double) * n);
fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nOut);
```
Note that above code is equivalent to the standard C way of allocating memory using `malloc`:
```C
double *in = (double*) malloc(sizeof(double) * n);
fftw_complex *out = (fftw_complex*) malloc(sizeof(fftw_complex) * nOut);
```
except that the FFTW routines ensure proper memory alignment for exploiting SIMD instructions of modern CPUs.

## Utility functions
In order to keep the examples short, a separate header file `util.h` is provided.
It contains methods to operate on one- and two-dimensional arrays (the latter stored in row-major order)
of real (`double`) and complex (`fftw_complex`) numbers.
The following operations are supported:
 * fill with random numbers between 0 and 1: e.g. `fill_random_1d_cplx`
 * element-wise check for approximate equality: e.g. `compare_1d_real`
 * write into a text file.

