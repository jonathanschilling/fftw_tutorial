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

int read_lcfs_header(char *filename, int *mpol, int *ntor, int *mnmax, int *nfp) {
    char *line = NULL;
    char *mpol_line = "# mpol\n";
    char *ntor_line = "# ntor\n";
    char *mnmax_line = "# mnmax\n";
    char *nfp_line = "# nfp\n";

    int found_mpol = 0;
    int found_ntor = 0;
    int found_mnmax = 0;
    int found_nfp = 0;

    size_t len = 0;
    ssize_t read;

    int status = 0;

    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("ERROR: could not open %st\n", filename);
        status = 1;
        return status;
    }

    int keep_going = 1;

    while (keep_going && (read = getline(&line, &len, fp)) != -1) {

        if (!found_mpol && strcmp(line, mpol_line) == 0) {
            read = getline(&line, &len, fp);
            if (read != -1) {
                (*mpol) = atoi(line);
                found_mpol = 1;
            } else {
                printf("failed to read line after '%s'\n", mpol_line);
                status = 1;
                keep_going = 0;
            }
        } else if (!found_ntor && strcmp(line, ntor_line) == 0) {
            read = getline(&line, &len, fp);
            if (read != -1) {
                (*ntor) = atoi(line);
                found_ntor = 1;
            } else {
                printf("failed to read line after '%s'\n", ntor_line);
                status = 1;
                keep_going = 0;
            }
        } else if (!found_mnmax && strcmp(line, mnmax_line) == 0) {
            read = getline(&line, &len, fp);
            if (read != -1) {
                (*mnmax) = atoi(line);
                found_mnmax = 1;
            } else {
                printf("failed to read line after '%s'\n", mnmax_line);
                status = 1;
                keep_going = 0;
            }
        } else if (!found_nfp && strcmp(line, nfp_line) == 0) {
            read = getline(&line, &len, fp);
            if (read != -1) {
                (*nfp) = atoi(line);
                found_nfp = 1;
            } else {
                printf("failed to read line after '%s'\n", nfp_line);
                status = 1;
                keep_going = 0;
            }
        }

        if (found_mpol && found_ntor && found_mnmax && found_nfp) {
            keep_going = 0;
        }
    }

    fclose(fp);

    if (line) {
        free(line);
    }

    return status;
}

int read_lcfs(char *filename, int mnmax, double *rmnc, double *zmns) {
    char *line = NULL;
    char *rmnc_line = "# rmnc\n";
    char *zmns_line = "# zmns\n";

    int found_rmnc = 0;
    int found_zmns = 0;

    size_t len = 0;
    ssize_t read;

    int status = 0;

    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("ERROR: could not open %s\n", filename);
        status = 1;
        return status;
    }

    int keep_going = 1;

    while (keep_going && (read = getline(&line, &len, fp)) != -1) {
        if (!found_rmnc && strcmp(line, rmnc_line) == 0) {
            for (int i=0; keep_going && i<mnmax; ++i) {
                read = getline(&line, &len, fp);
                if (read != -1) {
                    rmnc[i] = atof(line);
                } else {
                    printf("failed to read line %d after '%s'\n", i, rmnc_line);
                    status = 1;
                    keep_going = 0;
                }
            }
            found_rmnc = 1;
        } else if (!found_zmns && strcmp(line, zmns_line) == 0) {
            for (int i=0; keep_going && i<mnmax; ++i) {
                read = getline(&line, &len, fp);
                if (read != -1) {
                    zmns[i] = atof(line);
                } else {
                    printf("failed to read line %d after '%s'\n", i, zmns_line);
                    status = 1;
                    keep_going = 0;
                }
            }
            found_zmns = 1;
        }

        if (found_rmnc && found_zmns) {
            keep_going = 0;
        }
    }

    fclose(fp);

    if (line) {
        free(line);
    }

    return status;
}

#endif // UTIL_H
