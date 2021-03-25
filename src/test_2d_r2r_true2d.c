#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

int test_2d_r2r_true2d(int n0, int n1) {
    int idx_k, idx_j;

    double *in = fftw_alloc_real(n0 * n1);
    double *ref_out = fftw_alloc_real(n0 * n1);

    double *in1 = fftw_alloc_real(n0 * n1);
    double *in2 = fftw_alloc_real(n0 * n1);
    double *out1 = fftw_alloc_real(n0 * n1);
    double *out2 = fftw_alloc_real(n0 * n1);

    fftw_plan p1 = fftw_plan_r2r_2d(n0, n1, in1, out1, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
    fftw_plan p2 = fftw_plan_r2r_2d(n0, n1, in2, out2, FFTW_RODFT01, FFTW_RODFT01, FFTW_ESTIMATE);

    // random input
    fill_random_2d_real(n0, n1, in);

    // copy relevant input into in1 and in2
    int idx_j_2;
    for (int j0=0; j0<n0; ++j0) {
        for (int j1=0; j1<n1; ++j1) {
            idx_j = j0 * n1 + j1;

            // make sure that all entries are zero
            in1[idx_j] = 0.0;
            in2[idx_j] = 0.0;

            // input data for REDFT01
            double factor = 1.0;
            if (j0 > 0) {
                factor *= 0.5;
            }
            if (j1 > 0) {
                factor *= 0.5;
            }

            in1[idx_j] = factor * in[idx_j];

            // input data for RODFT01
            // need to shift indices by -1 due to (j+1) in definition
            if (j0 > 0 && j1 > 0) {
                int my_j0 = j0-1;
                int my_j1 = j1-1;

                idx_j_2 = my_j0 * n1 + my_j1;

                factor = 1.0;
                if (my_j0 < n0-1) {
                    factor *= 0.5;
                }
                if (my_j1 < n1-1) {
                    factor *= 0.5;
                }

                in2[idx_j_2] = factor * in[idx_j];
            }
        }
    }

    // manually compute 2D IDCT for reference
    double basis;
    for (int k0 = 0; k0 < n0; ++k0) {
        for (int k1 = 0; k1 < n1; ++k1) {
            idx_k = k0 * n1 + k1;

            ref_out[idx_k] = 0.0;

            for (int j0 = 0; j0 < n0; ++j0) {
                for (int j1 = 0; j1 < n1; ++j1) {
                    idx_j = j0 * n1 + j1;

                    basis = cos(M_PI * (  (k0 + 0.5) * j0 / ((double) n0)
                                        + (k1 + 0.5) * j1 / ((double) n1) ) );

                    ref_out[idx_k] += in[idx_j] * basis;
                }
            }
        }
    }

    fftw_execute(p1);
    fftw_execute(p2);

    // build combined output into out1
    for (int k0 = 0; k0 < n0; ++k0) {
        for (int k1 = 0; k1 < n1; ++k1) {
            idx_k = k0 * n1 + k1;

            // subtract because cos(u+v) = cos(u)*cos(v) - sin(u)*sin(v)
            out1[idx_k] -= out2[idx_k];
        }
    }

    // compare outputs
    double eps = 1.0e-12;
    int status = compare_2d_real(n0, n1, ref_out, out1, eps);

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_free(in);
    fftw_free(in1);
    fftw_free(in2);
    fftw_free(out1);
    fftw_free(out2);
    fftw_free(ref_out);

    return status;
}

int main(int argc, char **argv) {
    int status = 0;
    status += test_2d_r2r_true2d(4, 4);
    status += test_2d_r2r_true2d(4, 5);
    status += test_2d_r2r_true2d(5, 4);
    status += test_2d_r2r_true2d(5, 5);
    return status;
}
