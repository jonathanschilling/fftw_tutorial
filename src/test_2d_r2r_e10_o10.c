#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

int test_2d_r2r_e10_o10(int n0, int n1) {

    double *in = fftw_alloc_real(n0 * n1);
    double *out = fftw_alloc_real(n0 * n1);
    double *ref_out = fftw_alloc_real(n0 * n1);

    fftw_plan p = fftw_plan_r2r_2d(n0, n1, in, out, FFTW_REDFT10, FFTW_RODFT10, FFTW_ESTIMATE);

    // random input
    fill_random_2d_real(n0, n1, in);

    // manually compute DFT for reference
    int idx_k, idx_j;
    double basis_0, basis_1;
    for (int k0 = 0; k0 < n0; ++k0) {
        for (int k1 = 0; k1 < n1; ++k1) {
            idx_k = k0 * n1 + k1;

            ref_out[idx_k] = 0.0;

            for (int j0 = 0; j0 < n0; ++j0) {
                for (int j1 = 0; j1 < n1; ++j1) {
                    idx_j = j0 * n1 + j1;

                    // REDFT10 in first dimension
                    basis_0 = 2.0 * cos(M_PI * (j0 + 0.5) * k0 / ((double) n0));

                    // RODFT10 in second dimension
                    basis_1 = 2.0 * sin(M_PI * (j1 + 0.5) * (k1 + 1.0) / ((double) n1));

                    ref_out[idx_k] += in[idx_j] * basis_0 * basis_1;
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
    status += test_2d_r2r_e10_o10(4, 4);
    status += test_2d_r2r_e10_o10(4, 5);
    status += test_2d_r2r_e10_o10(5, 4);
    status += test_2d_r2r_e10_o10(5, 5);
    return status;
}
