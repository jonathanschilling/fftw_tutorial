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














}

int main(int argc, char** argv) {
    app_flux_surface(30, 36);

    return 0;
}
