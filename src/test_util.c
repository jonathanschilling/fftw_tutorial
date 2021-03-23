#include <complex.h>
#include <fftw3.h>

#include "util.h"

int main(int argc, char** argv) {

    int rows = 3;
    int cols = 4;

    fftw_complex *arr = fftw_alloc_complex(rows*cols);

    for (int i=0; i<rows; ++i) {
        for (int j=0; j<cols; ++j) {
            arr[i*cols+j] = i+I*j;
        }
    }

    dump_2d_cplx("test.dat", rows, cols, arr);

    // load in python as follows:
    // > import numpy as np
    // > d = np.loadtxt("test.dat", dtype=np.complex128)

    return 0;
}
