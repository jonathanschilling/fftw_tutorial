#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <complex.h>
#include <fftw3.h>

void fill_random_1d_real(int n, double *arr) {
    for (int i = 0; i < n; ++i) {
        arr[i] = rand() / ((double) RAND_MAX);
    }
}

void fill_random_1d_cplx(int n, fftw_complex *arr) {
    double real, imag;
    for (int i = 0; i < n; ++i) {
        real = rand() / ((double) RAND_MAX);
        imag = rand() / ((double) RAND_MAX);
        arr[i] = real + I * imag;
    }
}

void compare_1d_real(int n, double *ref, double *arr, double eps) {
    int status = 0;
    double delta;

    for (int i = 0; i < n; ++i) {
        printf("compare %d (len = %d)\n", i, n);

        delta = arr[i] - ref[i];

        if (delta < eps) {
            printf("  real ok (delta = %g)\n", delta);
        } else {
            printf("  real: expected %g, got %g (delta = %g)\n", ref[i], arr[i],
                    delta);
            status = 1;
        }
    }

    if (status == 0) {
        printf("=> all ok\n");
    } else {
        printf("=> errors\n");
    }
}

void compare_1d_cplx(int n, fftw_complex *ref, fftw_complex *arr, double eps) {
    int status = 0;
    double delta_real, delta_imag;

    for (int i = 0; i < n; ++i) {
        printf("compare %d (len = %d)\n", i, n);

        delta_real = creal(arr[i]) - creal(ref[i]);
        delta_imag = cimag(arr[i]) - cimag(ref[i]);

        if (delta_real < eps) {
            printf("  real ok (delta = %g)\n", delta_real);
        } else {
            printf("  real: expected %g, got %g (delta = %g)\n", creal(ref[i]),
                    creal(arr[i]), delta_real);
            status = 1;
        }

        if (delta_imag < eps) {
            printf("  imag ok (delta = %g)\n", delta_imag);
        } else {
            printf("  imag: expected %g, got %g (delta = %g)\n", cimag(ref[i]),
                    cimag(arr[i]), delta_imag);
            status = 1;
        }
    }

    if (status == 0) {
        printf("=> all ok\n");
    } else {
        printf("=> errors\n");
    }
}
//
//int compare_2d_cplx(int n_1, int n_2, fftw_complex *out, fftw_complex *ref_out,
//        double eps) {
//    int row_offset_k, index_k;
//    int status = 0;
//    double dR, dI;
//
//    row_offset_k = 0;
//    for (int k_1 = 0; k_1 < n_1; ++k_1) {
//        for (int k_2 = 0; k_2 < n_2; ++k_2) {
//            index_k = row_offset_k + k_2;
//
//            printf("compare ref_out[%d][%d] ...\n", k_1, k_2);
//
//            dR = creal(out[index_k]) - creal(ref_out[index_k]);
//            if (fabs(dR) > eps) {
//                printf("  real: expected %g, got %g (difference %g)\n",
//                        creal(out[index_k]), creal(ref_out[index_k]), dR);
//                status = 1;
//            } else {
//                printf("  real ok (%g)\n", dR);
//            }
//
//            dI = cimag(out[index_k]) - cimag(ref_out[index_k]);
//            if (fabs(dI) > eps) {
//                printf("  imag: expected %g, got %g (difference %g)\n",
//                        cimag(out[index_k]), cimag(ref_out[index_k]), dI);
//                status = 1;
//            } else {
//                printf("  imag ok (%g)\n", dI);
//            }
//
//        }
//        row_offset_k += n_2;
//    }
//    return status;
//}
//
//int compare_2d_real(int n_1, int n_2, double *out, double *ref_out, double eps) {
//    int row_offset_k, index_k;
//    int status = 0;
//    double d;
//
//    row_offset_k = 0;
//    for (int k_1 = 0; k_1 < n_1; ++k_1) {
//        for (int k_2 = 0; k_2 < n_2; ++k_2) {
//            index_k = row_offset_k + k_2;
//
//            printf("compare ref_out[%d][%d] ...\n", k_1, k_2);
//
//            d = out[index_k] - ref_out[index_k];
//            if (fabs(d) > eps) {
//                printf("  expected %g, got %g (difference %g)\n", out[index_k],
//                        ref_out[index_k], d);
//                status = 1;
//            } else {
//                printf("  ok (%g)\n", d);
//            }
//        }
//        row_offset_k += n_2;
//    }
//    return status;
//}
//
//void read_lcfs(int *mpol, int *ntor, int *mnmax, int *nfp, double *rmnc,
//        double *zmns) {
//    char *line = NULL;
//    char *comment = "#";
//    size_t comment_length = strlen(comment);
//    size_t len = 0;
//    ssize_t read;
//
//    FILE *fp = fopen("lcfs.dat", "r");
//    if (fp == NULL) {
//        printf("ERROR: could not open lcfs.dat\n");
//        return;
//    }
//
//    int i = 0;
//    while ((read = getline(&line, &len, fp)) != -1) {
//        // skip lines that start with '#'
//        if (strncmp(line, comment, comment_length) == 0) {
//            continue;
//        }
//
//        //printf("line %3d (len %3zu): %s", i, read, line);
//
//        // "parse" file contents
//        if (i == 0) {
//            (*mpol) = atoi(line);
//        } else if (i == 1) {
//            (*ntor) = atoi(line);
//        } else if (i == 2) {
//            (*mnmax) = atoi(line);
//        } else if (i == 3) {
//            (*nfp) = atoi(line);
//        } else if (i > 3 && i <= (*mnmax) + 3) {
//            rmnc[i - 4] = atof(line);
//        } else {
//            zmns[i - (*mnmax) - 4] = atof(line);
//        }
//
//        i++;
//    }
//
//    fclose(fp);
//
//    if (line) {
//        free(line);
//    }
//}
//
//void write_lcfs_realspace(char *filename, int nu, int nv, double *val) {
//
//    FILE *fp = fopen(filename, "w");
//
//    int index;
//    for (int i = 0; i < nv; ++i) {
//        for (int j = 0; j < nu; ++j) {
//            index = i * nu + j;
//
//            fprintf(fp, "%+6.5f ", val[index]);
//        }
//        fprintf(fp, "\n");
//    }
//
//    fclose(fp);
//}

#endif // UTIL_H
