#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

void fill_random_2d_real(int rows, int cols, double *arr) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            arr[i * cols + j] = rand() / ((double) RAND_MAX);
        }
    }
}

void fill_random_2d_cplx(int rows, int cols, fftw_complex *arr) {
    double real, imag;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            real = rand() / ((double) RAND_MAX);
            imag = rand() / ((double) RAND_MAX);
            arr[i * cols + j] = real + I * imag;
        }
    }
}

int compare_1d_real(int n, double *ref, double *arr, double eps) {
    int status = 0;
    double delta;

    for (int i = 0; i < n; ++i) {
        printf("compare %d (len = %d)\n", i, n);

        delta = arr[i] - ref[i];

        if (fabs(delta) < eps) {
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

    return status;
}

int compare_1d_cplx(int n, fftw_complex *ref, fftw_complex *arr, double eps) {
    int status = 0;
    double delta_real, delta_imag;

    for (int i = 0; i < n; ++i) {
        printf("compare %d (len = %d)\n", i, n);

        delta_real = creal(arr[i]) - creal(ref[i]);
        delta_imag = cimag(arr[i]) - cimag(ref[i]);

        if (fabs(delta_real) < eps) {
            printf("  real ok (delta = %g)\n", delta_real);
        } else {
            printf("  real: expected %g, got %g (delta = %g)\n", creal(ref[i]),
                    creal(arr[i]), delta_real);
            status = 1;
        }

        if (fabs(delta_imag) < eps) {
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

    return status;
}

int compare_2d_real(int rows, int cols, double *ref, double *arr, double eps) {
    int status = 0, index;
    double delta;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("compare (%d,%d) (shape = %dx%d)\n", i, j, rows, cols);

            index = i * cols + j;

            delta = arr[index] - ref[index];

            if (fabs(delta) < eps) {
                printf("  real ok (delta = %g)\n", delta);
            } else {
                printf("  real: expected %g, got %g (delta = %g)\n", ref[index],
                        arr[index], delta);
                status = 1;
            }
        }
    }

    if (status == 0) {
        printf("=> all ok\n");
    } else {
        printf("=> errors\n");
    }

    return status;
}

int compare_2d_cplx(int rows, int cols, fftw_complex *ref, fftw_complex *arr,
        double eps) {
    int status = 0, index;
    double delta_real, delta_imag;

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            printf("compare (%d,%d) (shape = %dx%d)\n", i, j, rows, cols);

            index = i * cols + j;

            delta_real = creal(arr[index]) - creal(ref[index]);
            delta_imag = cimag(arr[index]) - cimag(ref[index]);

            if (fabs(delta_real) < eps) {
                printf("  real ok (delta = %g)\n", delta_real);
            } else {
                printf("  real: expected %g, got %g (delta = %g)\n",
                        creal(ref[index]), creal(arr[index]), delta_real);
                status = 1;
            }

            if (fabs(delta_imag) < eps) {
                printf("  imag ok (delta = %g)\n", delta_imag);
            } else {
                printf("  imag: expected %g, got %g (delta = %g)\n",
                        cimag(ref[index]), cimag(arr[index]), delta_imag);
                status = 1;
            }
        }
    }

    if (status == 0) {
        printf("=> all ok\n");
    } else {
        printf("=> errors\n");
    }

    return status;
}

void dump_1d_real(char* filename, int n, double *arr) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("error: could not open file %s for writing\n", filename);
        return;
    }

    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%+6.5f\n", arr[i]);
    }

    fclose(fp);
}

void dump_1d_cplx(char* filename, int n, fftw_complex *arr) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("error: could not open file %s for writing\n", filename);
        return;
    }

    for (int i = 0; i < n; ++i) {
        fprintf(fp, "%+6.5f %+6.5f\n", creal(arr[i]), cimag(arr[i]));
    }

    fclose(fp);
}

void dump_2d_real(char* filename, int rows, int cols, double *arr) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("error: could not open file %s for writing\n", filename);
        return;
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            fprintf(fp, "%+6.5f ", arr[i * cols + j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

// using numpy format: matrix of (<real>+<cplx>j)
void dump_2d_cplx(char* filename, int rows, int cols, fftw_complex *arr) {
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("error: could not open file %s for writing\n", filename);
        return;
    }

    int index;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            index = i * cols + j;

            if (cimag(arr[index]) < 0) {
                fprintf(fp, "(%+6.5f%+6.5fj) ", creal(arr[index]), cimag(arr[index]));
            } else {
                fprintf(fp, "(%+6.5f+%6.5fj) ", creal(arr[index]), cimag(arr[index]));
            }
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

// manual implementation of a general DFT (shift and non-linear phase)
// formula from https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Generalized_DFT_(shifted_and_non-linear_phase) (accessed 2021-03-23)
// a: shift of input
// b: shift of output
void dft_1d_cplx(int n, fftw_complex *in, fftw_complex *out, double a, double b) {
    double phi;
    for (int k = 0; k < n; ++k) {
        out[k] = 0.0;
        for (int j = 0; j < n; ++j) {
            phi = -2.0 * M_PI * (k + b) * (j + a) / ((double) n);
            out[k] += in[j] * cexp(I * phi);
        }
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
