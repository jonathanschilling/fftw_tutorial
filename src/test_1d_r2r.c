#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

void test_1d_r2r() {

    int n = 5;

    double *in = fftw_alloc_real(n);
    double *out = fftw_alloc_real(n);

    fftw_plan p = fftw_plan_r2r_1d(n, in, out, FFTW_REDFT00, FFTW_ESTIMATE);

    fill_random_1d_real(n, in);

    dump_1d_real("test_1d_r2r_in.dat", n, in);

    fftw_execute(p);

    dump_1d_real("test_1d_r2r_out.dat", n, out);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}


int main(int argc, char** argv) {
    test_1d_r2r();

    return 0;
}
