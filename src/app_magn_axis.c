#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

void app_magn_axis(int n_zeta) {

    // number of Fourier coefficients
    int ntor = 12;

    // cos-parity Fourier coefficients for magnetic axis
    double R_ax_cos[13] = { 5.63, 0.391, 0.0123, 1.21e-3, 4.89e-6, -5.12e-5,
            -6.57e-5, 2.27e-6, -9.28e-5, -5.32e-7, 6.67e-5, 5.72e-5, 2.38e-5 };

    // sin-parity Fourier coefficients for magnetic axis
    double R_ax_sin[13] = { 0.0, 0.0727, 6.34e-3, 5.84e-3, 9.77e-4, 5.32e-5,
            8.48e-5, 5.57e-5, 5.56e-5, 5.53e-6, 7.74e-7, 1.03e-5, 8.75e-6 };

    int nCplx = n_zeta / 2 + 1;
    if (nCplx < ntor + 1) {
        printf("error: number of grid points too low.\n");
        return;
    }

    fftw_complex *in = fftw_alloc_complex(nCplx);
    double *out = fftw_alloc_real(n_zeta);

    fftw_plan p = fftw_plan_dft_c2r_1d(n_zeta, in, out, FFTW_ESTIMATE);

    // copy over available Fourier coefficients
    for (int n = 0; n <= ntor; ++n) {
        in[n] = R_ax_cos[n] - I * R_ax_sin[n];
    }

    // zero out remaining input Fourier coefficients
    for (int n = ntor + 1; n < nCplx; ++n) {
        in[n] = 0.0;
    }

    fftw_execute(p);

    // dump magnetic axis R geometry
    dump_1d_real("axis_R.dat", n_zeta, out);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

int main(int argc, char** argv) {
    app_magn_axis(36);

    return 0;
}
