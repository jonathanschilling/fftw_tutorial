#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

int test_1d_rodft00(int n) {

    int N = 2 * (n + 1);

    double *in = fftw_alloc_real(n);
    fftw_complex *in_logical = fftw_alloc_complex(N);
    double *out = fftw_alloc_real(n);
    fftw_complex *out_logical = fftw_alloc_complex(N);

    fftw_plan p = fftw_plan_r2r_1d(n, in, out, FFTW_RODFT00, FFTW_ESTIMATE);

    fill_random_1d_real(n, in);

    // the first half of the array is identical
    // except the first 0 (and the rest of the array being shifted)
    in_logical[0] = 0.0;
    for (int i = 0; i < n; ++i) {
        in_logical[i + 1] = in[i];
    }

    // second half is filled according to odd symmetry around (n+1)
    in_logical[n + 1] = 0.0;
    for (int i = 0; i < n; ++i) {
        in_logical[n + 2 + i] = -in[n - 1 - i];
    }

    dump_1d_real("test_1d_rodft00_in.dat", n, in);
    dump_1d_cplx("test_1d_rodft00_in_logical.dat", N, in_logical);

    // manual implementation of logically-equivalent DFT
    double a = 0.0;
    double b = 0.0;
    dft_1d_cplx(N, in_logical, out_logical, a, b);

    fftw_execute(p);

    dump_1d_real("test_1d_rodft00_out.dat", n, out);
    dump_1d_cplx("test_1d_rodft00_out_logical.dat", N, out_logical);

    // check output of logically equivalent DFT
    double eps = 1e-12;
    int status = 0;

    // 1. logically equivalent output should be purely imaginary-valued
    for (int i = 0; i < N; ++i) {
        if (fabs(creal(out_logical[i])) > eps) {
            printf("error: real of [%d] is %g\n", i, creal(out_logical[i]));
            status = 1;
        } else {
            printf("real of [%d] is %g\n", i, creal(out_logical[i]));
        }
    }

    // 2. odd symmetry about 0 implies that first entry should be zero
    if (fabs(cimag(out_logical[0])) > eps) {
        printf("error: delta of [%d] is %g\n", 0, -cimag(out_logical[0]));
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", 0, -cimag(out_logical[0]));
    }

    // 3. next n values should have the output of RODFT00 in the negative imaginary
    // with one index offset to account for the first zero in the input
    double delta;
    for (int i = 0; i < n; ++i) {
        delta = -cimag(out_logical[i + 1]) - out[i];
        if (fabs(delta) > eps) {
            printf("error: delta of [%d] is %g\n", i + 1, delta);
            status = 1;
        } else {
            printf("match of [%d] (delta=%g)\n", i + 1, delta);
        }
    }

    // 4. odd symmetry of output values around (n+1)
    if (fabs(cimag(out_logical[n + 1])) > eps) {
        printf("error: delta of [%d] is %g\n", n + 1,
                -cimag(out_logical[n + 1]));
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", n + 1, -cimag(out_logical[n + 1]));
    }
    for (int i = 0; i < n; ++i) {
        delta = -cimag(out_logical[n + 2 + i]) - (-out[n - 1 - i]);
        if (fabs(delta) > eps) {
            printf("error: delta of [%d] is %g\n", n + 2 + i, delta);
            status = 1;
        } else {
            printf("match of [%d] (delta=%g)\n", n + 2 + i, delta);
        }
    }

    if (status == 0) {
        printf("=> all ok\n");
    } else {
        printf("=> errors\n");
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(in_logical);
    fftw_free(out);
    fftw_free(out_logical);

    return status;
}

int main(int argc, char** argv) {
    int status = 0;
    status += test_1d_rodft00(3);
    status += test_1d_rodft00(4);
    return status;
}
