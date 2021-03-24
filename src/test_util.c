#include <complex.h>
#include <fftw3.h>

#include "util.h"

int main(int argc, char** argv) {

    int rows = 3;
    int cols = 4;

    fftw_complex *arr = fftw_alloc_complex(rows*cols);

    for (int i=0; i<rows; ++i) {
        for (int j=0; j<cols; ++j) {
            arr[i*cols+j] = i+I*(j-2);
        }
    }

    dump_2d_cplx("test_util.dat", rows, cols, arr);

    // load in python as follows:
    // > import numpy as np
    // > d = np.loadtxt("test_util.dat", dtype=np.complex128)
    // d should then have the following content:
    // array([[0.-2.j, 0.-1.j, 0.+0.j, 0.+1.j],
    //        [1.-2.j, 1.-1.j, 1.+0.j, 1.+1.j],
    //        [2.-2.j, 2.-1.j, 2.+0.j, 2.+1.j]])

    fftw_free(arr);

    return 0;
}
