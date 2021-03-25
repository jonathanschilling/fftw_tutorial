#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

void app_magn_axis_stellsym(int n_zeta) {

    // number of Fourier coefficients
    int ntor = 12;

    // cos-parity Fourier coefficients for R of magnetic axis
    double R_ax_cos[13] = { 5.63, 0.391, 0.0123, 1.21e-3, 4.89e-6, -5.12e-5,
            -6.57e-5, 2.27e-6, -9.28e-5, -5.32e-7, 6.67e-5, 5.72e-5, 2.38e-5 };

    int nCplx = n_zeta / 2;
    if (nCplx < ntor + 1) {
        printf("error: number of grid points too low.\n");
        return;
    }

    double *in = fftw_alloc_real(nCplx);
    double *out = fftw_alloc_real(n_zeta);

    fftw_plan p = fftw_plan_r2r_1d(nCplx, in, out, FFTW_REDFT01, FFTW_ESTIMATE);

    // copy over available Fourier coefficients
    in[0] = R_ax_cos[0];
    for (int n = 1; n <= ntor; ++n) {
        in[n] = 0.5 * R_ax_cos[n];
    }

    // zero out remaining input Fourier coefficients
    for (int n = ntor + 1; n < nCplx; ++n) {
        in[n] = 0.0;
    }

    fftw_execute(p);

    // fill second half of output array by mirroring about last entry
    for (int n = 0; n < nCplx; ++n) {
        out[n_zeta - 1 - n] = out[n];
    }

    // dump magnetic axis R geometry
    dump_1d_real("axis_R_halfGrid.dat", n_zeta, out);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
}

int main(int argc, char** argv) {
    app_magn_axis_stellsym(36);

    return 0;
}
