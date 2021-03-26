#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

void app_flux_surface(int n_theta, int n_zeta) {

    // read flux surface Fourier coefficients from file

    char *filename = "lcfs.dat";

    // number of poloidal modes
    int mpol;

    // number of toroidal modes
    int ntor;

    // total number of Fourier coefficients for R and L
    int mnmax;

    // number of field periods (5 for W7-X)
    int nfp;

    int status = 0;

    status = read_lcfs_header(filename, &mpol, &ntor, &mnmax, &nfp);
    if (status != 0) {
        printf("error in reading LCFS header\n");
        return;
    }

    printf(" mpol = %d\n", mpol);
    printf(" ntor = %d\n", ntor);
    printf("mnmax = %d\n", mnmax);
    printf("  nfp = %d\n", nfp);

    double *rmnc = fftw_alloc_real(mnmax);
    double *zmns = fftw_alloc_real(mnmax);

    status = read_lcfs(filename, mnmax, rmnc, zmns);
    if (status != 0) {
        printf("error in reading LCFS data\n");
        return;
    }

    // number of Fourier modes that can be represented with n_theta and n_zeta
    // without violating Nyquist criterion
    int nyq_pol = n_theta / 2 + 1;
    int nyq_tor = n_zeta / 2 + 1;

    if (nyq_pol < mpol) {
        printf("n_theta too small\n");
        return;
    }

    if (nyq_tor < ntor) {
        printf("n_zeta too small\n");
        return;
    }

    fftw_complex *in_R = fftw_alloc_complex(n_zeta * nyq_pol);
    fftw_complex *in_Z = fftw_alloc_complex(n_zeta * nyq_pol);
    double *out_R = fftw_alloc_real(n_zeta * n_theta);
    double *out_Z = fftw_alloc_real(n_zeta * n_theta);

    fftw_plan p_R = fftw_plan_dft_c2r_2d(n_zeta, n_theta, in_R, out_R, FFTW_ESTIMATE);
    fftw_plan p_Z = fftw_plan_dft_c2r_2d(n_zeta, n_theta, in_Z, out_Z, FFTW_ESTIMATE);

    fill_zero_2d_cplx(n_zeta, nyq_pol, in_R);
    fill_zero_2d_cplx(n_zeta, nyq_pol, in_Z);

    // copy LCFS Fourier coefficients into input array
    int idx_vmec = 0, idx_in;
    int m = 0;
    for (int n = 0; n <= ntor; ++n) {

        // compute target index for input array
        if (n <= 0) {
            idx_in = -n * nyq_pol + m;
        } else {
            idx_in = (n_zeta - n) * nyq_pol + m;
        }

        // VMEC: for m=0, only positive ntor --> summed up within VMEC definition with negative-n --> 2*0.5 = 1
        // also for sin-parity quantities: m=0 would imply (m theta -n zeta) == (-n zeta), but VMEC uses (n zeta) for m=0
        // --> does not matter for cos-parity since cos(-n zeta) = cos(n zeta), but sin(-n zeta) = -sin(n zeta)
        // and hence the additional -1 cancels out by compensating the -1 from the complex multiply in the DFT definition
        in_R[idx_in] = rmnc[idx_vmec];
        in_Z[idx_in] = I * zmns[idx_vmec];

        idx_vmec++;
    }

    for (m = 1; m < mpol; ++m) {
        for (int n = -ntor; n <= ntor; ++n) {

            // compute target index for input array
            if (n <= 0) {
                idx_in = -n * nyq_pol + m;
            } else {
                idx_in = (n_zeta - n) * nyq_pol + m;
            }

            // scale coefficients by 0.5 since FFTW implies
            // that coeffs for m<0 were present in the logically equivalent DFT input
            // which is not the case for VMEC output
            in_R[idx_in] = 0.5 * rmnc[idx_vmec];
            in_Z[idx_in] = -0.5 * I * zmns[idx_vmec];

            idx_vmec++;
        }
    }

    fftw_execute(p_R);
    fftw_execute(p_Z);

    dump_2d_real("lcfs_R.dat", n_zeta, n_theta, out_R);
    dump_2d_real("lcfs_Z.dat", n_zeta, n_theta, out_Z);

    fftw_destroy_plan(p_R);
    fftw_destroy_plan(p_Z);
    fftw_free(in_R);
    fftw_free(in_Z);
    fftw_free(out_R);
    fftw_free(out_Z);
}

int main(int argc, char** argv) {
    app_flux_surface(30, 36);

    return 0;
}
