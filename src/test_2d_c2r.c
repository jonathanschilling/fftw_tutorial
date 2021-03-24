#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

int test_2d_c2r(int n0, int n1) {
    int n1_cplx = n1 / 2 + 1;
    printf("n1_cplx = %d\n", n1_cplx);

    fftw_complex *in = fftw_alloc_complex(n0 * n1_cplx);
    double *out = fftw_alloc_real(n0 * n1);
    double *ref_out = fftw_alloc_real(n0 * n1);

    fftw_plan p = fftw_plan_dft_c2r_2d(n0, n1, in, out, FFTW_ESTIMATE);

    // random input
    fill_random_2d_cplx(n0, n1_cplx, in);

    // manually compute DFT for reference
    int idx_k, idx_j;
    double phi, real;
    for (int k0 = 0; k0 < n0; ++k0) {
        for (int k1 = 0; k1 < n1; ++k1) {
            idx_k = k0 * n1 + k1;

            ref_out[idx_k] = 0.0;

            for (int j0 = 0; j0 < n0; ++j0) {
                for (int j1 = 0; j1 < n1_cplx; ++j1) {
                    idx_j = j0 * n1_cplx + j1;

                    phi = 2.0 * M_PI * (  k0 * j0 / ((double) n0)
                                        + k1 * j1 / ((double) n1) );

                    // output is purely real,
                    // so compute only real part
                    real = creal(in[idx_j]) * cos(phi) - cimag(in[idx_j]) * sin(phi);

                    ref_out[idx_k] += real;

                    // add symmetric entries twice
                    // n1/2+n1%2 is n1/2 if n1 is even
                    // and it is n1/2+1 if n1 is odd
                    if (j1 > 0 && j1 < n1 / 2 + n1 % 2) {
                        ref_out[idx_k] += real;
                    }
                }
            }
        }
    }

    fftw_execute(p);

    // compare outputs
    double eps = 1.0e-12;
    int status = compare_2d_real(n0, n1, ref_out, out, eps);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_free(ref_out);

    return status;
}

int main(int argc, char **argv) {
    int status = 0;
    status += test_2d_c2r(4, 4);
    status += test_2d_c2r(4, 5);
    status += test_2d_c2r(5, 4);
    status += test_2d_c2r(5, 5);
    return status;
}
