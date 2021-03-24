#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

int test_2d_c2c(int n0, int n1) {
    int dir = FFTW_FORWARD;

    fftw_complex *in = fftw_alloc_complex(n0 * n1);
    fftw_complex *out = fftw_alloc_complex(n0 * n1);
    fftw_complex *ref_out = fftw_alloc_complex(n0 * n1);

    fftw_plan p = fftw_plan_dft_2d(n0, n1, in, out, dir, FFTW_ESTIMATE);

    // random input
    fill_random_2d_cplx(n0, n1, in);

    // manually compute DFT for reference
    int idx_k, idx_j;
    double phi;
    for (int k0 = 0; k0 < n0; ++k0) {
        for (int k1 = 0; k1 < n1; ++k1) {
            idx_k = k0 * n1 + k1;

            ref_out[idx_k] = 0.0;

            for (int j0 = 0; j0 < n0; ++j0) {
                for (int j1 = 0; j1 < n1; ++j1) {
                    idx_j = j0 * n1 + j1;

                    phi = dir * 2.0 * M_PI * (
                              k0 * j0 / ((double) n0)
                            + k1 * j1 / ((double) n1) );

                    ref_out[idx_k] += in[idx_j] * cexp(I * phi);
                }
            }
        }
    }

    fftw_execute(p);

    // compare outputs
    double eps = 1.0e-12;
    int status = compare_2d_cplx(n0, n1, ref_out, out, eps);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_free(ref_out);

    return status;
}

int main(int argc, char **argv) {
    int status = 0;
    status += test_2d_c2c(4, 4);
    status += test_2d_c2c(4, 5);
    status += test_2d_c2c(5, 4);
    status += test_2d_c2c(5, 5);
    return status;
}

