#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

void test_1d_c2r(int n) {
    double real, phi;

    int nCplx = n / 2 + 1;

    printf("n=%d nCplx=%d\n", n, nCplx);

    fftw_complex *in = fftw_alloc_complex(nCplx);
    double *ref_out = fftw_alloc_real(n);
    double *fftw_out = fftw_alloc_real(n);

    fftw_plan p = fftw_plan_dft_c2r_1d(n, in, fftw_out, FFTW_ESTIMATE);

    // fill the input array with random data
    fill_random_1d_cplx(nCplx, in);

    // compute the reference output
    for (int k = 0; k < n; ++k) {

        // start with DC component, which is purely real due to Hermitian symmetry
        ref_out[k] = creal(in[0]);

        if (n % 2 == 0) {
            // even n
            for (int j = 1; j < nCplx - 1; ++j) {
                phi = 2.0 * M_PI * j * k / ((double) n);

                real = creal(in[j]) * cos(phi) - cimag(in[j]) * sin(phi);
                ref_out[k] += 2.0 * real;
            }

            // Nyquist element is purely real as well
            phi = 2.0 * M_PI * (nCplx - 1) * k / ((double) n);
            ref_out[k] += creal(in[nCplx - 1]) * cos(phi);

        } else {
            // odd n
            for (int j = 1; j < nCplx; ++j) {
                phi = 2.0 * M_PI * j * k / ((double) n);

                real = creal(in[j]) * cos(phi) - cimag(in[j]) * sin(phi);
                ref_out[k] += 2.0 * real;
            }
        }
    }

    // compute the DFT of in using FFTW
    fftw_execute(p);

    // compare reference output with FFTW output
    double eps = 1e-12;
    compare_1d_real(n, ref_out, fftw_out, eps);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(ref_out);
    fftw_free(fftw_out);
}

int main(int argc, char** argv) {
    test_1d_c2r(32);
    test_1d_c2r(33);
    return 0;
}