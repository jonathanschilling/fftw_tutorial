# FFTW Tutorial

This is a basic C project (Makefile, but also for the Eclipse IDE) I use for exploring FFTW 3.3.9.

One- and two-dimensional discrete Fourier transforms (DFTs) of random data are computed using both FFTW and straight-forward naive algorithms
in order to illustrate explicitly what kind of symmetries and scaling properties FFTW implies in its inputs and outputs.

## Table of Contents:
1. [One-Dimensional Examples](https://github.com/jonathanschilling/fftw_tutorial#one-dimensional-examples)
   1. [1D complex-to-complex](https://github.com/jonathanschilling/fftw_tutorial#1d-complex-to-complex)
   2. [1D complex-to-real and real-to-complex](https://github.com/jonathanschilling/fftw_tutorial#1d-complex-to-real-and-real-to-complex)
   3. [1D real-to-real](https://github.com/jonathanschilling/fftw_tutorial#1d-real-to-real)
      1. [REDFT00 (DCT-I)](https://github.com/jonathanschilling/fftw_tutorial#redft00-dct-i)
      2. [REDFT10 (DCT-II)](https://github.com/jonathanschilling/fftw_tutorial#redft10-dct-ii)
      3. [REDFT01 (DCT-III)](https://github.com/jonathanschilling/fftw_tutorial#redft01-dct-iii)
      4. [REDFT11 (DCT-IV)](https://github.com/jonathanschilling/fftw_tutorial#redft11-dct-iv)
      5. [RODFT00 (DST-I)](https://github.com/jonathanschilling/fftw_tutorial#rodft00-dst-i)
      6. [RODFT10 (DST-II)](https://github.com/jonathanschilling/fftw_tutorial#rodft10-dst-ii)
      7. [RODFT01 (DST-III)](https://github.com/jonathanschilling/fftw_tutorial#rodft01-dst-iii)
      8. [RODFT11 (DST-IV)](https://github.com/jonathanschilling/fftw_tutorial#rodft11-dst-iv)
2. [Two-Dimensional Examples](https://github.com/jonathanschilling/fftw_tutorial#two-dimensional-examples)
   1. [2D complex-to-complex](https://github.com/jonathanschilling/fftw_tutorial#2d-complex-to-complex)
   2. [2D complex-to-real](https://github.com/jonathanschilling/fftw_tutorial#2d-complex-to-real)
   3. [2D real-to-complex](https://github.com/jonathanschilling/fftw_tutorial#2d-real-to-complex)
   4. [2D real-to-real](https://github.com/jonathanschilling/fftw_tutorial#2d-real-to-real)
      1. [2D REDFT10+RODFT10](https://github.com/jonathanschilling/fftw_tutorial#2d-redft10rodft10)
      2. [True 2D IDCT using FFTW](https://github.com/jonathanschilling/fftw_tutorial#true-2d-idct-using-fftw)
3. [Real-world applications](https://github.com/jonathanschilling/fftw_tutorial#real-world-applications)
   1. [Geometry of the Magnetic Axis in a Stellarator](https://github.com/jonathanschilling/fftw_tutorial#geometry-of-the-magnetic-axis-in-a-stellarator)
   2. [Geometry of a Flux Surface in a Stellarator](https://github.com/jonathanschilling/fftw_tutorial#geometry-of-a-flux-surface-in-a-stellarator)
4. [Allocation of arrays](https://github.com/jonathanschilling/fftw_tutorial#allocation-of-arrays)
5. [Utility functions](https://github.com/jonathanschilling/fftw_tutorial#utility-functions)

## One-Dimensional Examples
This tutorial starts by computing one-dimensional (1D) DFTs of random input data.

### 1D complex-to-complex
The first example is basically a self-contained version
of the [corresponding example in the FFTW manual](http://fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html#Complex-One_002dDimensional-DFTs). 

We want to compute the complex-valued one-dimensional DFT here, which is specified in
[section 4.8.1 of the FFTW reference manual](http://fftw.org/fftw3_doc/The-1d-Discrete-Fourier-Transform-_0028DFT_0029.html#The-1d-Discrete-Fourier-Transform-_0028DFT_0029).

![complex-valued DFT](eqn/complex_dft.png)

The sign in the exponent of the basis function specifies the direction in which the Fourier transform is to be computed:
`-1` indicates a "forward" transform and `+1` indicates a backward transform.
These values are available via the `FFTW_FORWARD` and `FFTW_BACKWARD` preprocessor macros.

In order to compute the DFT, complex-valued products of the following form need to be evaluated:

![complex product](eqn/c2c_product.png)

Eulers formula comes in handy now (where *i* is the imaginary unit with *i*^2=-1):

![Eulers formula](eqn/euler.png)

The angle argument `phi` can be identified in above formulas for the DFT:

![phi in DFT](eqn/phi.png)

Now the complex-valued product can be computed using only real-valued variables:

![complex product using real numbers](eqn/c2c_product_real.png)

FFTW implements all this math internally and the explicit formulation was only used to build a basis for the computations to follow.
Below is the example code showing how to compute the 1d `c2c` DFT using both FFTW and a manual implementation.
The size of the DFT is specified via the variable `n` and the direction (forward or backward) is specified via the variable `dir`.
Complex-valued arrays `in`, `ref_out` and `fftw_out` are allocated to hold the input (*X_k*) and outputs (*Y_k*) of the DFT.
A plan for the corresponding DFT is created using `fftw_plan_dft_1d`.
Only after this, the input data is written into the `in` array.
The reference output is computed now (before calling `fftw_execute`), since in some later examples (most notably multi-dimensional `c2r` transforms),
FFTW overwrites the input data and for good practise, we keep this in mind already now.
Note that this is an out-of-place transform, since `in` and `fftw_out` are allocated to be separate arrays.
Next the FFTW transform can be executed via `fftw_execute`.
This fills the corresponding output array `fftw_out`, which is subsequently compared against the reference output in `ref_out`.
A conservative tolerance of `1e-12` is specified to make the example work also for weird input data (as generated by the PRNG).
The actual deviations are usually much smaller and can be observed in the screen output from `compare_1d_cplx`.
Finally, the `fftw_plan` is destroyed and the memory is released using `fftw_free` (which has to be used if the array was allocated using `fftw_alloc_*`).

```C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

int test_1d_c2c(int n) {

    int dir = FFTW_BACKWARD;
    double real, imag, phi;

    fftw_complex *in = fftw_alloc_complex(n);
    fftw_complex *ref_out = fftw_alloc_complex(n);
    fftw_complex *fftw_out = fftw_alloc_complex(n);

    fftw_plan p = fftw_plan_dft_1d(n, in, fftw_out, dir, FFTW_ESTIMATE);

    // fill the input array with random data
    fill_random_1d_cplx(n, in);

    // compute the reference output
    for (int k = 0; k < n; ++k) {
        ref_out[k] = 0.0;
        for (int j = 0; j < n; ++j) {
            phi = dir * 2.0 * M_PI * j * k / ((double) n);

            real = creal(in[j]) * cos(phi) - cimag(in[j]) * sin(phi);
            imag = creal(in[j]) * sin(phi) + cimag(in[j]) * cos(phi);
            ref_out[k] += real + I * imag;
        }
    }

    // compute the DFT of in using FFTW
    fftw_execute(p);

    // compare reference output with FFTW output
    double eps = 1e-12;
    int status = compare_1d_cplx(n, ref_out, fftw_out, eps);

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(ref_out);
    fftw_free(fftw_out);

    return status;
}

int main(int argc, char** argv) {
    int status = 0;
    status += test_1d_c2c(32);
    status += test_1d_c2c(33);
    return status;
}
```

The code is available in the file [`src/test_1d_c2c.c`](src/test_1d_c2c.c).

### 1D complex-to-real and real-to-complex
The next two examples deal with DFTs of purely real data (`r2c`) and DFTs which *produce* purely real data (`c2r`).
These are covered in the [official FFTW tutorial](http://fftw.org/fftw3_doc/One_002dDimensional-DFTs-of-Real-Data.html#One_002dDimensional-DFTs-of-Real-Data)
as well as in the [FFTW reference manual](http://fftw.org/fftw3_doc/The-1d-Real_002ddata-DFT.html#The-1d-Real_002ddata-DFT).

In case either the input array or the output array are constrained to be purely real, the corresponding complex-valued output or input array
features Hermitian symmetry (where the `n`-periodicity has been included as well):

![Hermitian symmetry](eqn/hermitian.png)

For the case of a DFT of real-valued *X_j* and complex-valued *Y_k* with Hermitian symmetry,
the Fourier sum is written out excplicitly as follows:

![Hermitian symmetry in Fourier sum](eqn/hermitian_sum.png)

The figure below illustrates the structure of the complex-valued Fourier space arrays
occuring in the DFT for both even-valued (`n=6`) and odd-valued (`n=7`) sizes of the DFT.

![array structure](img/array_structures.png)

The size required to contain all information required for the transform from or to a real-valued array
is contained in the first `n/2+1` (division by 2 rounded down) entries of the complex array, indicated by the red bars in above figure.
For both even and odd values of `n`, Hermitian symmetry implies *Y_0* = *Y\*_0* and thus *Y_0* is always real.
In the case of even `n`, we can intuitively observe that *Y_m* = *Y\*_-m* where *m*=`n/2` (the element at the Nyquist frequency,
and thus also the last element of the complex-valued Fourier-space array) is also purely real.

The DFT formulation includes all elements of the Fourier-space array from 0 to `n-1`.
Now that only parts of these coefficients are taken into account, they have to be weighted appropriately
to recover the results that would have been obtained by using the full array without taking advantage of its symmetry properties.
For odd `n`, all components of the Fourier-space array except the DC element at *k*=0 have to be weighted with a factor of 2.
For even `n`, all components of the Fourier-space array except the DC element at *k*=0 and the Nyquist element at *k*=`n/2`
have to be weighted with a factor of 2.
The elements that need to be weighted by a factor of 2 are highlighted by solid blue lines in above illustration.
The redundant elements that are not explicitly needed are indicated by dashed blue lines.

The backward transformation from complex-valued Fourier space
to real space is demonstrated in [`src/test_1d_c2r.c`](src/test_1d_c2r.c).
The relevant portion of the source code is here:

```C
int nCplx = n / 2 + 1;
for (int k = 0; k < n; ++k) {

    // start with DC component, which is purely real due to Hermitian symmetry
    ref_out[k] = creal(in[0]);

    int loopEnd = nCplx;

    // special case for even n
    if (n % 2 == 0) {
        // Nyquist element is purely real as well
        phi = 2.0 * M_PI * (nCplx - 1) * k / ((double) n);
        ref_out[k] += creal(in[nCplx - 1]) * cos(phi);

        loopEnd = nCplx-1;
    }

    // middle elements are handled the same for even and odd n
    for (int j = 1; j < loopEnd; ++j) {
        phi = 2.0 * M_PI * j * k / ((double) n);

        real = creal(in[j]) * cos(phi) - cimag(in[j]) * sin(phi);
        ref_out[k] += 2.0 * real;
    }
}
```

Note that the integer division used to compute `nCplx` (correctly) rounds down.
The rest is a relatively straight-forward implementation of above verbose algorithm.
The DC component is always taken to be real.
Depending on whether `n` is even or odd, the number of elements to take into account with both real and imaginary component (`loopEnd`) is adjusted.
The (purely real) Nyquist element at `n/2` is added separately if `n` is even.
All other elements are weighted by a factor of 2 and only the real part
of the complex product of input Fourier coefficient and complex-valued basis function is actually computed.

The forward transform from real space to Fourier space is comparably simple to implement:

```C
for (int k = 0; k < nCplx; ++k) {

    // DC component is always real
    ref_out[k] = in[0];

    for (int j = 1; j < n; ++j) {
        phi = -2.0 * M_PI * j * k / ((double) n);

        real = in[j] * cos(phi);
        imag = in[j] * sin(phi);
        ref_out[k] += real + I * imag;
    }
}
```

Note that in this case, only the non-redundant part of the complex-values Fourier coefficients need to be computed
from a real-valued input. The separate handling of the DC component is not strictly necessary, since `cos(0)=1` and `sin(0)=0`
and thus the DC component would get no imaginary contribution.
The full example is available in the file [`src/test_1d_r2c.c`](src/test_1d_r2c.c).

It becomes evident in above examples that the sign of a `c2r` DFT is always `FFTW_BACKWARD` and 
the sign of a `r2c` DFT is always `FFTW_FORWARD`.

### 1D real-to-real

Certain symmetries can be assumed for a given real input array of which a DFT is to be computed
that lead to the output array being purely real as well. This is another gain of a factor of 2 in speed
and memory usage over `r2c`/`c2r` transforms.
Depending on even (odd) parity of the input array, the transform outputs have even (odd) parity.
They are therefore called Discrete Cosine Transform (DCT) and Discrete Sine Transfrom (DST), respectively.

The logical size of the corresponding DFT is denoted as *N*.
The actual array size given to FFTW is denoted by `n`.
For the DCT/DST types implemented in FFTW, *N* is always even.
Note that this does not pose any restrictions on the input array sizes `n`.
One can think of the logical DFT input array as one that FFTW 'sees' internally and computes a regular DFT of.
The resulting output array is purely real and features the named symmetry properties,
since the logical input array was 'constructed' from the given input array to have the desired symmetry properties.

In below example plots used to illustrate these symmetry properties/assumptions,
random Fourier coefficients have been sampled and transformed back to real-space using the named inverse transforms.
This allows to 'evaluate' the input data also in the samples given in the (small) input arrays.
In these plots, red dashed vertical lines indicate even symmetry (*f(x)=f(-x)*) about the indicated position
and blue dashed vertical lines indicate odd symmetry (*f(x)=-f(-x)*) about the indicated position.
The grey-shaded area in the background indicated the range of samples that are included in the input array.
The x axis labels denote the indices in the actual input array given to FFTW (only valid within the range of the grey boxes).

The nomenclature works as follows:
The first letter is **R** to indicate real-valued data.
The second letter distinguished between **E** for even-parity data and **O** for odd-parity data.
The following **DFT** is for discrete Fourier transform (who guessed...).
The next two digits indicate wheter (**1**) or not (**0**) the input (first digit) or the output (second digit) data is shifted by half a sample.
Think of this in terms of parity: whether the symmetry axis is located at a sample (no shifting necessary) or between two samples (shifting necessary).
The shifting becomes necessary when formulating the symmetry properties over sampled data that has integer indices
vs. symmetry axis that are possibly located at half-integer locations.
The parity and symmetry properties of the output array are those of the input array for the inverse transform.

For all transforms, a periodicity of *N* is assumed for the *logical* input array as *X_j = X_{N+j}* where *X* is the input data array.

Here is a quick overview table to indicate the assumed symmetries in the input array for the following types of `r2r` DFTs:

| type    | actual `r2r` input | logically-equivalent DFT input |
| :------:| :---------: | :-------------------: |
| REDFT00 | `a b c d e` | `a b c d  e  d  c  b` |
| REDFT10 | `a b c d  ` | `a b c d  d  c  b  a` |
| REDFT01 | `a b c d  ` | `a b c d  0 -d -c -b` |
| REDFT11 | `a b c d  ` | `a b c d -d -c -b -a` |
| RODFT00 | `a b c    ` | `0 a b c  0 -c -b -a` |
| RODFT10 | `a b c d  ` | `a b c d -d -c -b -a` |
| RODFT01 | `a b c d  ` | `a b c d  c  b  a  0` |
| RODFT11 | `a b c d  ` | `a b c d  d  c  b  a` |

For the DCTs, please also consider https://en.wikipedia.org/wiki/Discrete_cosine_transform#Formal_definition
and in particular https://upload.wikimedia.org/wikipedia/commons/a/ae/DCT-symmetries.svg .

For the DSTs, please also consider https://en.wikipedia.org/wiki/Discrete_sine_transform#Definition
and in particular https://upload.wikimedia.org/wikipedia/commons/3/31/DST-symmetries.svg .

#### REDFT00 (DCT-I)

In case of the real-valued even-parity DFT with no shifts in either input or output array (REDFT00),
also called the DCT-I, the corresponding logical DFT size is given by *N* = 2(`n`-1), corresponding to `n` = *N*/2+1.

The formal definition of the REDFT00 is given below:

![REDFT00 formula](eqn/redft00.png)

The inverse of this transform is REDFT00 itself.
The input array is assumed to have even symmetry around *j=0* and even symmetry also around *j=n−1*.

![REDFT00](img/redft00.png)

In above figure, the lowercase letters *a* to *e* refer to the input data *abcde* for the size-5 REDFT00,
which is logically equivalent to a size-8 DFT with real-valued input data *abcdedcb*.

In order to demonstrate the use of this method,
the logically equivalent DFT input is filled appropriately and its output is checked against `REDFT00`.
In the following code, `in` is the input array (size `n`) given to `REDFT00`
and `in_logical` is the (complex-valued) input array (size *N*) handed to a
[generic 1D DFT](https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Generalized_DFT_(shifted_and_non-linear_phase)).
Similarly, `out` is the output array (size `n`) from `REDFT00`
and `out_logical` is the output array (size *N*) from a generic 1D DFT.

Here is how the symmetric input is generated:

```C
// the first half of the array is identical
for (int i=0; i<n; ++i) {
    in_logical[i] = in[i];
}

// second half is filled according to even symmetry around n-1
for (int i=0; i<n-1; ++i) {
    in_logical[n+i] = in[n-2-i];
}
```

The checks are a little bit more involved.
The logically equivalent DFT output should be purely real-valued:

```C
for (int i = 0; i < N; ++i) {
    if (fabs(cimag(out_logical[i])) > eps) {
        printf("error: imag of [%d] is %g\n", i, cimag(out_logical[i]));
        status = 1;
    } else {
        printf("imag of [%d] is %g\n", i, cimag(out_logical[i]));
    }
}
```

The first `n` values should be identical between `REDFT00` and the generalized DFT:

```C
double delta;
for (int i = 0; i < n; ++i) {
    delta = creal(out_logical[i]) - out[i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", i, delta);
    }
}
```

The remaining values should have even symmetry around `n-1`:

```C
for (int i = 0; i < n - 2; ++i) {
    delta = creal(out_logical[n + i]) - out[n - 2 - i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", n+i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", n+i, delta);
    }
}
```

The full example can be found in [`src/test_1d_redft00.c`](src/test_1d_redft00.c).

#### REDFT10 (DCT-II)

In case of the real-valued even-parity DFT with shifted input data (REDFT10),
also called the DCT-II, the corresponding logical DFT size is given by *N* = 2`n`, corresponding to `n` = *N*/2.
This function is commonly known as "the" DCT.

The formal definition of the REDFT10 is given below:

![REDFT10 formula](eqn/redft10.png)

The inverse of this transform is REDFT01.
The input array is assumed to have even symmetry around *j=-0.5* and even symmetry also around *j=n−0.5*.

![REDFT10](img/redft10.png)

In above figure, the lowercase letters *a* to *e* refer to the input data *abcd* for the size-4 REDFT10,
which is logically equivalent to a size-8 DFT with real-valued input data *abcddcba*.

In order to demonstrate the use of this method,
the logically equivalent DFT input is filled appropriately and its output is checked against `REDFT10`.
In the following code, `in` is the input array (size `n`) given to `REDFT10`
and `in_logical` is the (complex-valued) input array (size *N*) handed to a
[generic 1D DFT](https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Generalized_DFT_(shifted_and_non-linear_phase)).
Similarly, `out` is the output array (size `n`) from `REDFT10`
and `out_logical` is the output array (size *N*) from a generic 1D DFT.

Here is how the symmetric input is generated:

```C
// the first half of the array is identical
for (int i = 0; i < n; ++i) {
    in_logical[i] = in[i];
}

// second half is filled according to even symmetry around n-0.5
for (int i = 0; i < n; ++i) {
    in_logical[n + i] = in[n - 1 - i];
}
```

The checks are a little bit more involved.
The logically equivalent DFT output should be purely real-valued:

```C
for (int i = 0; i < N; ++i) {
    if (fabs(cimag(out_logical[i])) > eps) {
        printf("error: imag of [%d] is %g\n", i, cimag(out_logical[i]));
        status = 1;
    } else {
        printf("imag of [%d] is %g\n", i, cimag(out_logical[i]));
    }
}
```

The first `n` values should be identical between `REDFT10` and the generalized DFT:

```C
for (int i = 0; i < n; ++i) {
    delta = creal(out_logical[i]) - out[i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", i, delta);
    }
}
```

Odd symmetry around `n` implies that the value at `n` is zero:

```C
if (fabs(creal(out_logical[n])) > eps) {
    printf("error: delta of [%d] is %g\n", n, creal(out_logical[n]));
    status = 1;
} else {
    printf("match of [%d] (delta=%g)\n", n, creal(out_logical[n]));
}
```

The remaining values should have odd symmetry around `n`:

```C
for (int i = 1; i < n; ++i) {
    delta = creal(out_logical[n + i]) - (-out[n - i]);
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", n + i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", n + i, delta);
    }
}
```

The full example can be found in [`src/test_1d_redft10.c`](src/test_1d_redft10.c).

#### REDFT01 (DCT-III)

In case of the real-valued even-parity DFT with shifted output data (REDFT01),
also called the DCT-III, the corresponding logical DFT size is given by *N* = 2`n`, corresponding to `n` = *N*/2.
This function is commonly known as "the" inverse DCT (IDCT).

The formal definition of the REDFT01 is given below:

![REDFT01 formula](eqn/redft01.png)

The inverse of this transform is REDFT10.
The input array is assumed to have even symmetry around *j=0* and odd symmetry around *j=n*.

![REDFT01](img/redft01.png)

In order to demonstrate the use of this method,
the logically equivalent DFT input is filled appropriately and its output is checked against `REDFT01`.
In the following code, `in` is the input array (size `n`) given to `REDFT01`
and `in_logical` is the (complex-valued) input array (size *N*) handed to a
[generic 1D DFT](https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Generalized_DFT_(shifted_and_non-linear_phase)).
Similarly, `out` is the output array (size `n`) from `REDFT01`
and `out_logical` is the output array (size *N*) from a generic 1D DFT.

Here is how the symmetric input is generated:

```C
// the first half of the array is identical
for (int i = 0; i < n; ++i) {
    in_logical[i] = in[i];
}

// second half is filled according to odd symmetry around n
for (int i = 1; i < n; ++i) {
    in_logical[n + i] = -in[n - i];
}
```

The checks are a little bit more involved.
The logically equivalent DFT output should be purely real-valued:

```C
for (int i = 0; i < N; ++i) {
    if (fabs(cimag(out_logical[i])) > eps) {
        printf("error: imag of [%d] is %g\n", i, cimag(out_logical[i]));
        status = 1;
    } else {
        printf("imag of [%d] is %g\n", i, cimag(out_logical[i]));
    }
}
```

The first `n` values should be identical between `REDFT01` and the generalized DFT:

```C
for (int i = 0; i < n; ++i) {
    delta = creal(out_logical[i]) - out[i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", i, delta);
    }
}
```

The remaining values should have even symmetry around `n-0.5`:

```C
for (int i = 0; i < n; ++i) {
    delta = creal(out_logical[n + i]) - out[n - 1 - i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", n + i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", n + i, delta);
    }
}
```

The full example can be found in [`src/test_1d_redft01.c`](src/test_1d_redft01.c).

#### REDFT11 (DCT-IV)

In case of the real-valued even-parity DFT with both shifted input and output data (REDFT11),
also called the DCT-IV, the corresponding logical DFT size is given by *N* = 2`n`, corresponding to `n` = *N*/2.

The formal definition of the REDFT11 is given below:

![REDFT11 formula](eqn/redft11.png)

The inverse of this transform is REDFT11 itself.
The input array is assumed to have even symmetry around *j=-0.5* and odd symmetry around *j=n-0.5*.

![REDFT11](img/redft11.png)

In order to demonstrate the use of this method,
the logically equivalent DFT input is filled appropriately and its output is checked against `REDFT11`.
In the following code, `in` is the input array (size `n`) given to `REDFT11`
and `in_logical` is the (complex-valued) input array (size *N*) handed to a
[generic 1D DFT](https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Generalized_DFT_(shifted_and_non-linear_phase)).
Similarly, `out` is the output array (size `n`) from `REDFT11`
and `out_logical` is the output array (size *N*) from a generic 1D DFT.

Here is how the symmetric input is generated:

```C
// the first half of the array is identical
for (int i = 0; i < n; ++i) {
    in_logical[i] = in[i];
}

// second half is filled according to odd symmetry around n-0.5
for (int i = 0; i < n; ++i) {
    in_logical[n + i] = -in[n - 1 - i];
}
```

The checks are a little bit more involved.
The logically equivalent DFT output should be purely real-valued:

```C
for (int i = 0; i < N; ++i) {
    if (fabs(cimag(out_logical[i])) > eps) {
        printf("error: imag of [%d] is %g\n", i, cimag(out_logical[i]));
        status = 1;
    } else {
        printf("imag of [%d] is %g\n", i, cimag(out_logical[i]));
    }
}
```

The first `n` values should be identical between `REDFT11` and the generalized DFT:

```C
for (int i = 0; i < n; ++i) {
    delta = creal(out_logical[i]) - out[i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", i, delta);
    }
}
```

The remaining values should have odd symmetry around `n-0.5`:

```C
for (int i = 0; i < n; ++i) {
    delta = creal(out_logical[n + i]) - (-out[n - 1 - i]);
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", n + i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", n + i, delta);
    }
}
```

The full example can be found in [`src/test_1d_redft11.c`](src/test_1d_redft11.c).

#### RODFT00 (DST-I)

In case of the real-valued odd-parity DFT with no shifts in either input or output array (RODFT00),
also called the DST-I, the corresponding logical DFT size is given by *N* = 2(`n`+1), corresponding to `n` = *N*/2-1.
Note that the periodicity of *N* of the logical input array in combination with odd symmetry *X_j = -X_{n-j}*
leads to *X_0 = -X_0* which is equivalent to *X_0 = 0*.
This first always-zero element of the input array is not explicitly included in the input to FFTW
and the input array thus has a size of one less and the indices of the symmetry axis shift by 1.

The formal definition of the RODFT00 is given below:

![RODFT00 formula](eqn/rodft00.png)

The inverse of this transform is RODFT00 itself.
The input array is assumed to have odd symmetry around *j=-1* and odd symmetry also around *j=n*.

![RODFT00](img/rodft00.png)

In order to demonstrate the use of this method,
the logically equivalent DFT input is filled appropriately and its output is checked against `RODFT00`.
In the following code, `in` is the input array (size `n`) given to `RODFT00`
and `in_logical` is the (complex-valued) input array (size *N*) handed to a
[generic 1D DFT](https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Generalized_DFT_(shifted_and_non-linear_phase)).
Similarly, `out` is the output array (size `n`) from `RODFT00`
and `out_logical` is the output array (size *N*) from a generic 1D DFT.

Here is how the symmetric input is generated:

```C
// the first half of the array is identical
// except the first 0 (and the rest of the array being shifted)
in_logical[0] = 0.0;
for (int i = 0; i < n; ++i) {
    in_logical[i + 1] = in[i];
}

// second half is filled according to odd symmetry around (n+1)
in_logical[n + 1] = 0.0;
for (int i = 0; i < n; ++i) {
    in_logical[n + 2 + i] = -in[n - 1 - i];
}
```

The checks are a little bit more involved.
The logically equivalent DFT output should be purely imaginary-valued:

```C
for (int i = 0; i < N; ++i) {
    if (fabs(creal(out_logical[i])) > eps) {
        printf("error: real of [%d] is %g\n", i, creal(out_logical[i]));
        status = 1;
    } else {
        printf("real of [%d] is %g\n", i, creal(out_logical[i]));
    }
}
```

Odd symmetry about `0` implies that the imaginary part of the first output should be zero:

```C
if (fabs(cimag(out_logical[0])) > eps) {
    printf("error: delta of [%d] is %g\n", 0, -cimag(out_logical[0]));
    status = 1;
} else {
    printf("match of [%d] (delta=%g)\n", 0, -cimag(out_logical[0]));
}
```

The next `n` values should have the output of `RODFT00` in their negative imaginary parts
with one index offset to account for the first zero in the input:

```C
for (int i = 0; i < n; ++i) {
    delta = -cimag(out_logical[i + 1]) - out[i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", i + 1, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", i + 1, delta);
    }
}
```

Odd symmetry about `n` in the input to `RODFT00` implies
odd symmetry in the output of the logically-equivalent DFT about `n+1` due to the index shift by one:

```C
if (fabs(cimag(out_logical[n + 1])) > eps) {
    printf("error: delta of [%d] is %g\n", n + 1,
            -cimag(out_logical[n + 1]));
    status = 1;
} else {
    printf("match of [%d] (delta=%g)\n", n + 1, -cimag(out_logical[n + 1]));
}
```

The remaining values should have odd symmetry around `n+1`
(note again that the real output from `RODFT00` needs to be compared against
the negative imaginary parts in the output of the logically equivalent DFT shifted by one full index):

```C
for (int i = 0; i < n; ++i) {
    delta = -cimag(out_logical[n + 2 + i]) - (-out[n - 1 - i]);
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", n + 2 + i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", n + 2 + i, delta);
    }
}
```

The full example can be found in [`src/test_1d_rodft00.c`](src/test_1d_rodft00.c).

#### RODFT10 (DST-II)

In case of the real-valued odd-parity DFT with shifted input data (RODFT10),
also called the DST-II, the corresponding logical DFT size is given by *N* = 2`n`, corresponding to `n` = *N*/2.
This function is commonly known as "the" DST.

The formal definition of the RODFT10 is given below:

![RODFT10 formula](eqn/rodft10.png)

The inverse of this transform is RODFT01.
The input array is assumed to have odd symmetry around *j=-0.5* and odd symmetry also around *j=n-0.5*.

![RODFT10](img/rodft10.png)

In order to demonstrate the use of this method,
the logically equivalent DFT input is filled appropriately and its output is checked against `RODFT10`.
In the following code, `in` is the input array (size `n`) given to `RODFT10`
and `in_logical` is the (complex-valued) input array (size *N*) handed to a
[generic 1D DFT](https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Generalized_DFT_(shifted_and_non-linear_phase)).
Similarly, `out` is the output array (size `n`) from `RODFT10`
and `out_logical` is the output array (size *N*) from a generic 1D DFT.

Here is how the symmetric input is generated:

```C
// the first half of the array is identical
for (int i = 0; i < n; ++i) {
    in_logical[i] = in[i];
}

// second half is filled according to odd symmetry around (n-0.5)
for (int i = 0; i < n; ++i) {
    in_logical[n + i] = -in[n - 1 - i];
}
```

The checks are a little bit more involved.
The logically equivalent DFT output should be purely imaginary-valued:

```C
for (int i = 0; i < N; ++i) {
    if (fabs(creal(out_logical[i])) > eps) {
        printf("error: real of [%d] is %g\n", i, creal(out_logical[i]));
        status = 1;
    } else {
        printf("real of [%d] is %g\n", i, creal(out_logical[i]));
    }
}
```

The output of `RODFT10` is shifted by one index to the left with respect to the logically-equivalent DFT,
since the first input is constrained to be zero (resulting in odd parity about `-1` in the input to `RODFT10`).
Odd symmetry about `0` implies that the imaginary part of the first output should be zero:

```C
if (fabs(cimag(out_logical[0])) > eps) {
    printf("error: delta of [%d] is %g\n", 0, -cimag(out_logical[0]));
    status = 1;
} else {
    printf("match of [%d] (delta=%g)\n", 0, -cimag(out_logical[0]));
}
```

The next `n` values should have the output of `RODFT10` in their negative imaginary parts
with one index offset to account for the first zero in the input:

```C
for (int i = 0; i < n; ++i) {
    delta = -cimag(out_logical[i + 1]) - out[i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", i + 1, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", i + 1, delta);
    }
}
```

The remaining values have even symmetry around `n-1`
(note again that the real output from `RODFT10` needs to be compared against
the negative imaginary parts in the output of the logically equivalent DFT shifted by one full index):

```C
for (int i = 0; i < n - 1; ++i) {
    delta = -cimag(out_logical[n + 1 + i]) - out[n - 2 - i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", n + 1 + i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", n + 1 + i, delta);
    }
}
```

The full example can be found in [`src/test_1d_rodft10.c`](src/test_1d_rodft10.c).

#### RODFT01 (DST-III)

In case of the real-valued odd-parity DFT with shifted output data (RODFT01),
also called the DST-III, the corresponding logical DFT size is given by *N* = 2`n`, corresponding to `n` = *N*/2.
This function is commonly known as "the" inverse DST (IDST).

The formal definition of the RODFT01 is given below:

![RODFT01 formula](eqn/rodft01.png)

The inverse of this transform is RODFT10.
The input array is assumed to have odd symmetry around *j=-1* and even symmetry around *j=n-1*.

![RODFT01](img/rodft01.png)

In order to demonstrate the use of this method,
the logically equivalent DFT input is filled appropriately and its output is checked against `RODFT01`.
In the following code, `in` is the input array (size `n`) given to `RODFT01`
and `in_logical` is the (complex-valued) input array (size *N*) handed to a
[generic 1D DFT](https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Generalized_DFT_(shifted_and_non-linear_phase)).
Similarly, `out` is the output array (size `n`) from `RODFT01`
and `out_logical` is the output array (size *N*) from a generic 1D DFT.

The input of `RODFT01` is shifted by one index to the left with respect to the logically-equivalent DFT,
since the first input is constrained to be zero (resulting in odd parity about `-1` in the input to `RODFT01`).
Here is how the symmetric input is generated:

```C
// a zero in the first entry is needed to satisfy
// odd symmetry about -1 in the shifted input array
in_logical[0] = 0.0;

// the first half of the array is identical up to shift
for (int i = 0; i < n; ++i) {
    in_logical[i + 1] = in[i];
}

// second half is filled according to even symmetry around (n-1)
for (int i = 0; i < n - 1; ++i) {
    in_logical[n + 1 + i] = in[n - 2 - i];
}
```

The checks are a little bit more involved.
The logically equivalent DFT output should be purely imaginary-valued:

```C
for (int i = 0; i < N; ++i) {
    if (fabs(creal(out_logical[i])) > eps) {
        printf("error: real of [%d] is %g\n", i, creal(out_logical[i]));
        status = 1;
    } else {
        printf("real of [%d] is %g\n", i, creal(out_logical[i]));
    }
}
```

The first `n` values should have the output of `RODFT01` in their negative imaginary parts:

```C
for (int i = 0; i < n; ++i) {
    delta = -cimag(out_logical[i]) - out[i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", i, delta);
    }
}
```

The remaining values have odd symmetry around `n-0.5`
(note again that the real output from `RODFT01` needs to be compared against
the negative imaginary parts in the output of the logically equivalent DFT shifted by one full index):

```C
for (int i = 0; i < n; ++i) {
    delta = -cimag(out_logical[n + i]) - (-out[n - 1 - i]);
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", n + i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", n + i, delta);
    }
}
```

The full example can be found in [`src/test_1d_rodft01.c`](src/test_1d_rodft01.c).

#### RODFT11 (DST-IV)

In case of the real-valued odd-parity DFT with both shifted input and output data (RODFT11),
also called the DST-IV, the corresponding logical DFT size is given by *N* = 2`n`, corresponding to `n` = *N*/2.

The formal definition of the RODFT11 is given below:

![RODFT11 formula](eqn/rodft11.png)

The inverse of this transform is RODFT11 itself.
The input array is assumed to have odd symmetry around *j=-0.5* and even symmetry around *j=n-0.5*.

![RODFT11](img/rodft11.png)

In order to demonstrate the use of this method,
the logically equivalent DFT input is filled appropriately and its output is checked against `RODFT11`.
In the following code, `in` is the input array (size `n`) given to `RODFT11`
and `in_logical` is the (complex-valued) input array (size *N*) handed to a
[generic 1D DFT](https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Generalized_DFT_(shifted_and_non-linear_phase)).
Similarly, `out` is the output array (size `n`) from `RODFT11`
and `out_logical` is the output array (size *N*) from a generic 1D DFT.

Here is how the symmetric input is generated:

```C
// the first half of the array is identical up to shift
for (int i = 0; i < n; ++i) {
    in_logical[i] = in[i];
}

// second half is filled according to even symmetry around (n-0.5)
for (int i = 0; i < n; ++i) {
    in_logical[n + i] = in[n - 1 - i];
}
```

The checks are a little bit more involved.
The logically equivalent DFT output should be purely imaginary-valued:

```C
for (int i = 0; i < N; ++i) {
    if (fabs(creal(out_logical[i])) > eps) {
        printf("error: real of [%d] is %g\n", i, creal(out_logical[i]));
        status = 1;
    } else {
        printf("real of [%d] is %g\n", i, creal(out_logical[i]));
    }
}
```

The first `n` values should have the output of `RODFT11` in their negative imaginary parts:

```C
for (int i = 0; i < n; ++i) {
    delta = -cimag(out_logical[i]) - out[i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", i, delta);
    }
}
```

The remaining values have even symmetry around `n-0.5`
(note again that the real output from `RODFT11` needs to be compared against
the negative imaginary parts in the output of the logically equivalent DFT):

```C
for (int i = 0; i < n; ++i) {
    delta = -cimag(out_logical[n + i]) - out[n - 1 - i];
    if (fabs(delta) > eps) {
        printf("error: delta of [%d] is %g\n", n + i, delta);
        status = 1;
    } else {
        printf("match of [%d] (delta=%g)\n", n + i, delta);
    }
}
```

The full example can be found in [`src/test_1d_rodft11.c`](src/test_1d_rodft11.c).

## Two-Dimensional Examples
The next section deals with the computation of two-dimensional (2D) DFTs of random input data.

### 2D complex-to-complex

A first example for the computation of two-dimensional Fourier transforms using FFTW is given here:

```C
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <complex.h>
#include <fftw3.h>

#include "util.h"

int test_2d_c2c(int n0, int n1) {
    int dir = FFTW_FORWARD;

    fftw_complex *in = fftw_alloc_complex(n0 * n1);
    fftw_complex *out = fftw_alloc_complex(n0 * n1);
    fftw_complex *ref_out = fftw_alloc_complex(n0 * n1);

    fftw_plan p = fftw_plan_dft_2d(n0, n1, in, out, dir, FFTW_ESTIMATE);

    // random input
    fill_random_2d_cplx(n0, n1, in);

    // manually compute DFT for reference
    int idx_k, idx_j;
    double phi;
    for (int k0 = 0; k0 < n0; ++k0) {
        for (int k1 = 0; k1 < n1; ++k1) {
            idx_k = k0 * n1 + k1;

            ref_out[idx_k] = 0.0;

            for (int j0 = 0; j0 < n0; ++j0) {
                for (int j1 = 0; j1 < n1; ++j1) {
                    idx_j = j0 * n1 + j1;

                    phi = dir * 2.0 * M_PI * (
                              k0 * j0 / ((double) n0)
                            + k1 * j1 / ((double) n1) );

                    ref_out[idx_k] += in[idx_j] * cexp(I * phi);
                }
            }
        }
    }

    fftw_execute(p);

    // compare outputs
    double eps = 1.0e-12;
    int status = compare_2d_cplx(n0, n1, ref_out, out, eps);

    return status;
}

int main(int argc, char **argv) {
    int status = 0;
    status += test_2d_c2c(4, 4);
    status += test_2d_c2c(4, 5);
    status += test_2d_c2c(5, 4);
    status += test_2d_c2c(5, 5);
    return status;
}
```

As usual, arrays are allocated to the specified size and filled with random data.
The two-dimensional complex-valued DFT is computed manually for reference
and the outputs are compared against the result of FFTW.
The product of the one-dimensional DFTs along the two dimensions is computed here
via addition of the corresponding phases in the exponent of the complex-valued basis function.

The code can be found in [`src/test_2d_c2c.c`](src/test_2d_c2c.c).

### 2D complex-to-real

The two-dimensional `c2r` transform can make use of the Hermitian symmetry
([see above](https://github.com/jonathanschilling/fftw_tutorial#1d-complex-to-real-and-real-to-complex))
to reduce the computational work in the *last* dimension of the input data by about a factor of 2.
The input data thus has to have a shape of (`n0`x`n1_cplx`) where `n1_cplx = n1/2+1` (division by 2 rounded down)
and (`n0`x`n1`) is the logical size of the DFT and also the size of the real-valued output arrays.

The arrays are allocated as follows:

```C
fftw_complex *in = fftw_alloc_complex(n0 * n1_cplx);
double *out = fftw_alloc_real(n0 * n1);
double *ref_out = fftw_alloc_real(n0 * n1);
```

The manual transform that is computed for reference is done as follows:

```C
int idx_k, idx_j;
double phi, real;
for (int k0 = 0; k0 < n0; ++k0) {
    for (int k1 = 0; k1 < n1; ++k1) {
        idx_k = k0 * n1 + k1;

        ref_out[idx_k] = 0.0;

        for (int j0 = 0; j0 < n0; ++j0) {
            for (int j1 = 0; j1 < n1_cplx; ++j1) {
                idx_j = j0 * n1_cplx + j1;

                phi = 2.0 * M_PI * (  k0 * j0 / ((double) n0)
                                    + k1 * j1 / ((double) n1) );

                // output is purely real,
                // so compute only real part
                real = creal(in[idx_j]) * cos(phi) - cimag(in[idx_j]) * sin(phi);

                ref_out[idx_k] += real;

                // add symmetric entries twice
                // n1/2+n1%2 is n1/2 if n1 is even
                // and it is n1/2+1 if n1 is odd
                if (j1 > 0 && j1 < n1 / 2 + n1 % 2) {
                    ref_out[idx_k] += real;
                }
            }
        }
    }
}
```

The full example can be found in [`src/test_2d_c2r.c`](src/test_2d_c2r.c).

### 2D real-to-complex

The two-dimensional `r2c` transform can make use of the Hermitian symmetry
([see above](https://github.com/jonathanschilling/fftw_tutorial#1d-complex-to-real-and-real-to-complex))
to reduce the computational work in the *last* dimension of the output data by about a factor of 2.
The real-valued input data has to have a shape of (`n0`x`n1`), which is also the logical size of the DFT.
The output data then has a shape of (`n0`x`n1_cplx`) where `n1_cplx = n1/2+1` (division by 2 rounded down).

The arrays are allocated as follows:

```C
double *in = fftw_alloc_real(n0 * n1);
fftw_complex *out = fftw_alloc_complex(n0 * n1_cplx);
fftw_complex *ref_out = fftw_alloc_complex(n0 * n1_cplx);
```

The manual transform that is computed for reference is done as follows:

```C
int idx_k, idx_j;
double phi;
for (int k0 = 0; k0 < n0; ++k0) {
    for (int k1 = 0; k1 < n1_cplx; ++k1) {
        idx_k = k0 * n1_cplx + k1;

        ref_out[idx_k] = 0.0;

        for (int j0 = 0; j0 < n0; ++j0) {
            for (int j1 = 0; j1 < n1; ++j1) {
                idx_j = j0 * n1 + j1;

                phi = -2.0 * M_PI * (  k0 * j0 / ((double) n0)
                                     + k1 * j1 / ((double) n1) );

                ref_out[idx_k] += in[idx_j] * cexp(I*phi);
            }
        }
    }
}
```

The full example can be found in [`src/test_2d_r2c.c`](src/test_2d_r2c.c).

### 2D real-to-real

Two-dimensional real-to-real DFTs can be defined with any of the one-dimensional transforms
in each of the two dimensions. This implies that for eight one-dimensional transforms
there is a total of 64 possible combinations that exceed the scope of this tutorial
and also do not contribute much insight.

Note that FFTW computes the *separable product* of the transforms along each dimension
and not the true multi-dimensional DFT.
This can be understood intuitively by the addition formulas for sine and cosine:

![Addition formulas for sine and cosine](eqn/sin_cos_add.png)

FFTW's two-dimensional real-to-real DFTs can be used to compute the products on the right-hand side
of above equation and two such DFts have to be added together to arrive at the true
two-dimensional DFT with the Fourier kernel as on the left-hand side of above equation.

In the following, two examples are presented.
First, a mixed `r2r` DFT with `REDFT10` in the first dimension and `RODFT10` in the second dimension is computed.
Second, a true two-dimensional DFT in computed using two FFTW calls.

#### 2D REDFT10+RODFT10

An example of a two-dimensional transform is presented
with `REDFT10` in the first dimension and `RODFT10` in the second dimension.

The basis functions for the two dimensions are computed separately:

```C
// REDFT10 in first dimension
basis_0 = 2.0 * cos(M_PI * (j0 + 0.5) * k0 / ((double) n0));

// RODFT10 in second dimension
basis_1 = 2.0 * sin(M_PI * (j1 + 0.5) * (k1 + 1.0) / ((double) n1));
```

The contribution to the output array is computed as follows:

```C
ref_out[idx_k] += in[idx_j] * basis_0 * basis_1;
```

The full example can be found in [`src/test_2d_r2r_e10_o10.c`](src/test_2d_r2r_e10_o10.c).

#### True 2D IDCT using FFTW

A true two-dimensional inverse discrete cosine transform is computed using FFTW in this example:

![True two-dimensional IDCT](eqn/true_2d_idct.png)

FFTW can only compute the separable product of one-dimensional transforms for multi-dimensional transform.
Above kernel function thus has to be split up by application of the addition formula for sine and cosine:

![True two-dimensional IDCT kernel](eqn/true_2d_idct_kernel.png)

The left half of above equation is computed using a two-dimensional `REDFT01`
and the right half is computed using a two-dimensional `RODFT01`.
The outputs of these transforms are subtracted from each other to yield the desired result.

First consider the naive implementation in plain `C`.

The input and output arrays are allocated and filled with random data:

```C
double *in = fftw_alloc_real(n0 * n1);
double *ref_out = fftw_alloc_real(n0 * n1);

fill_random_2d_real(n0, n1, in);
```

The two-dimensional IDCT is computed from the input data:

```C
double basis;
for (int k0 = 0; k0 < n0; ++k0) {
    for (int k1 = 0; k1 < n1; ++k1) {
        idx_k = k0 * n1 + k1;

        ref_out[idx_k] = 0.0;

        for (int j0 = 0; j0 < n0; ++j0) {
            for (int j1 = 0; j1 < n1; ++j1) {
                idx_j = j0 * n1 + j1;

                basis = cos(M_PI * (  (k0 + 0.5) * j0 / ((double) n0)
                                    + (k1 + 0.5) * j1 / ((double) n1) ) );

                ref_out[idx_k] += in[idx_j] * basis;
            }
        }
    }
}
```

Next, the same transform is implemented using FFTW.

Separate input and output array are allocated in this educational example:

```C
double *in1 = fftw_alloc_real(n0 * n1);
double *in2 = fftw_alloc_real(n0 * n1);
double *out1 = fftw_alloc_real(n0 * n1);
double *out2 = fftw_alloc_real(n0 * n1);
```

The two two-dimensional transforms are planned:

```C
fftw_plan p1 = fftw_plan_r2r_2d(n0, n1, in1, out1, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
fftw_plan p2 = fftw_plan_r2r_2d(n0, n1, in2, out2, FFTW_RODFT01, FFTW_RODFT01, FFTW_ESTIMATE);
```

The next step is to fill the two input arrays `in1` and `in2` for FFTW
with appropriate data so that the results of FFTW are equivalent
to those of the naive implementation presented above.
This step is not entirely trivial and presented in detail next.

A two-dimensional loop over the input data `in` is set up
and all inputs to FFTW are initialized to be zero.
Then, the two input arrays are filled separately.

```C
for (int j0=0; j0<n0; ++j0) {
    for (int j1=0; j1<n1; ++j1) {
        idx_j = j0 * n1 + j1;

        // make sure that all entries are zero
        in1[idx_j] = 0.0;
        in2[idx_j] = 0.0;
        
        /* fill in1 */
        
        /* fill in2 */
    }
}
```

The input data for `REDFT01` has to be multiplied by `0.5`
for each dimension where the index is above `0` to cancel
the factor of `2` in the [definition of `REDFT01`](https://github.com/jonathanschilling/fftw_tutorial#redft01-dct-iii)
for direct correspondence with the desired true two-dimensional IDCT.
This is the `/* fill in1 */` step in above loop:

```C
double factor = 1.0;
if (j0 > 0) {
    factor *= 0.5;
}
if (j1 > 0) {
    factor *= 0.5;
}

in1[idx_j] = factor * in[idx_j];
```

The input data for `RODFT01` also needs to have entries
multiplied by `0.5` for each dimension in which the index is less than `n-1`
to cancel the factor of `2` in the [definition of `RODFT01`](https://github.com/jonathanschilling/fftw_tutorial#rodft01-dst-iii)
for direct correspondence with the desired true two-dimensional IDCT.
Additionally, the input data needs to be shifted left/down by one index
due to the `j+1` appearing in the definition of `RODFT01`.
This is the `/* fill in2 */` step in above loop:

```C
if (j0 > 0 && j1 > 0) {
    int my_j0 = j0-1;
    int my_j1 = j1-1;

    idx_j_2 = my_j0 * n1 + my_j1;

    factor = 1.0;
    if (my_j0 < n0-1) {
        factor *= 0.5;
    }
    if (my_j1 < n1-1) {
        factor *= 0.5;
    }

    in2[idx_j_2] = factor * in[idx_j];
}
```

In the desired truely-2D IDCT, there is no contribution from the second term
in the split-basis approach used to implement this transform using FFTW 
if either `j0==0` or `j1==0`, since the respecitve arguments to `sin` are zero in this case.
This motivates the check for `j0 > 0 && j1 > 0` above.
The indices `j0` and `j1` in the `REDFT01` input are transformed to the corresponding indices `my_j0` and `my_j1`
in the input to `RODFT01`. The corresponding linear index `idx_j_2` in `in2` is computed from `my_j0` and `my_j1`.
Factors of `0.5` are required for each dimension in which the index is less than `n-1`
to cancel the factor of `2` in the definition of `RODFT01`. 

The FFTW plans can be executed after computing the reference output
and it remains to combine the outputs by subtracting the contribution
from the two-dimensional `RODFT01` to arrive at the final result in `out1`:

```C
for (int k0 = 0; k0 < n0; ++k0) {
    for (int k1 = 0; k1 < n1; ++k1) {
        idx_k = k0 * n1 + k1;

        // subtract because cos(u+v) = cos(u)*cos(v) - sin(u)*sin(v)
        out1[idx_k] -= out2[idx_k];
    }
}
```

The full example can be found in [`src/test_2d_r2r_true2d.c`](src/test_2d_r2r_true2d.c).

## Real-world applications

In this section, FFTW is used to evaluate the geometry
of the magnetic axis and a flux surface of a stellarator ideal magnetohydrodynamic (MHD) plasma equilibrium
as computed by the [Variational Moments Equilibrium Code (VMEC)](https://doi.org/10.1016/0010-4655(86)90058-5).

Here is a plot of a few of the flux surfaces of the [Wendelstein 7-X](https://www.ipp.mpg.de/w7x) stellarator:

![flux surfaces and magnetic axis of W7-X](img/flux_surfaces_w7x_ref_163_1_10_50_98_cut_small.png)

Nested flux surfaces are shown in grey, red, green and blue.
The black line indicates the magnetic axis.

The geometry of the magnetic axis and that of a flux surface are real-valued.
Therefore, `c2r` DFTs are sufficient to compute the backward transforms to real space.

### Geometry of the Magnetic Axis in a Stellarator
The real-space geometry of the magnetic axis (a general closed curve) is given via a one-dimensional DFT
[evaluated using real-valued arithmetics](https://github.com/jonathanschilling/fftw_tutorial#1d-complex-to-complex):

![full Fourier series for magnetic axis](eqn/magn_axis.png)

where ζ is the toroidal angle per field period in radians.
For a five-fold symmetric stellarator like W7-X, ζ ranges from 0 to 5 * 2π around the whole machine.
Conversely, the first (unique) toroidal segment of the geometry is contained within a range of ζ from 0 to 2π.

In this example the geometry of the magnetic axis is evaluated at regular intervals in ζ:

![grid in zeta](eqn/grid_zeta.png)

where *l* ranges from `0` to `n_ζ-1`.

Inserting this into the Fourier series for, e.g., the *R* coordinate, leads to a DFT:

![magnetic axis R via DFT](eqn/magn_axis_dft.png)

Note that typically the number of Fourier coefficients included in computation
of the MHD equilibrium is (much) less than the number of grid points, i.e., *N* < `n_ζ`.
The physics happens in real space and certains aspects of a computation
might require a fine spatial resolution (i.e., many grid points in real space) to accurately
resolve, e.g., gradients while the overall complexity of the solution
(gouverned by the number of included Fourier coefficients) does not need to be extremely high.
FFTW only supports (logically) equally-sized input and output arrays and thus,
if less Fourier coefficients than desired grid points are to be used,
the Fourier coefficient array needs to be filled with zeros in the remaining area.

The axis geometry is defined as follows:

```C
// number of Fourier coefficients
int ntor = 12;

// cos-parity Fourier coefficients for magnetic axis
double R_ax_cos[13] = { 5.63, 0.391, 0.0123, 1.21e-3, 4.89e-6, -5.12e-5,
    -6.57e-5, 2.27e-6, -9.28e-5, -5.32e-7, 6.67e-5, 5.72e-5, 2.38e-5 };

// sin-parity Fourier coefficients for magnetic axis
double R_ax_sin[13] = { 0.0, 0.0727, 6.34e-3, 5.84e-3, 9.77e-4, 5.32e-5,
    8.48e-5, 5.57e-5, 5.56e-5, 5.53e-6, 7.74e-7, 1.03e-5, 8.75e-6 };
```

The number of grid points in the toroidal direction at which the magnetic axis geometry
is to be evaluated is specified via the `n_zeta` parameter.
It needs to be checked that `n_zeta` is large enough, i.e.,
above the Nyquist limit corresponding to the specified number of Fourier coefficients.
Then, the input and output arrays can be allocated:

```C
int nCplx = n_zeta/2+1;
if (nCplx<ntor+1) {
    printf("error: number of grid points too low.\n");
    return;
}

fftw_complex *in = fftw_alloc_complex(nCplx);
double *out = fftw_alloc_real(n_zeta);

fftw_plan p = fftw_plan_dft_c2r_1d(n_zeta, in, out, FFTW_ESTIMATE);
```

Now that the input array is allocated,
the available Fourier coefficients can be copied over
and (if present) the remainder of the input array is set to zero:

```C
// copy over available Fourier coefficients
for (int n=0; n<=ntor; ++n) {
    in[n] = R_ax_cos[n] - I * R_ax_sin[n];
}

// zero out remaining input Fourier coefficients
for (int n=ntor+1; n<nCplx; ++n) {
    in[n] = 0.0;
}
```

The full example can be found in [`src/app_magn_axis.c`](src/app_magn_axis.c).

Assuming [stellarator symmetry](https://doi.org/10.1016/S0167-2789(97)00216-9),
half of the Fourier coefficients can be omitted and the transform
reduces to the one-dimensional IDCT and IDST, respectively:

![stellarator-symmetric Fourier series for magnetic axis](eqn/magn_axis_stellsym.png)

The [FFTW documentation](http://fftw.org/fftw3_doc/Real-even_002fodd-DFTs-_0028cosine_002fsine-transforms_0029.html#DOCF4)
explicitly advises against using `R*DFT00` transforms for this purpose due to numerical stability issues
which are currently circumvented within FFTW by using a less efficient algorithm.
If the data is required anyway on the whole domain from 0 to 2π,
it is probably easiest and fastest (in terms of development work)
to simply use a `c2r` DFT and provide only real data in the input.

On the other hand, in case the evaluation locations are not required to be located on
"full" grid points in the toroidal direction, `REDFT01` and  `RODFT01` can be used.
This is the subject of the now-starting second half of this example.

The number of grid points as well as the output array size are given as follows:

```C
int nCplx = n_zeta / 2;
if (nCplx < ntor + 1) {
    printf("error: number of grid points too low.\n");
    return;
}
```

Only the cosine Fourier coefficients are copied over into the (real) input array:

```C
// copy over available Fourier coefficients
in[0] = R_ax_cos[0];
for (int n = 1; n <= ntor; ++n) {
    in[n] = R_ax_cos[n];
}

// zero out remaining input Fourier coefficients
for (int n = ntor + 1; n < nCplx; ++n) {
    in[n] = 0.0;
}
```

The second half of the output array (of size `n_zeta` for equivalent with the first half of this example)
is filled by applying the even symmetry of the output data implicitly used in `REDFT01`:

```C
for (int n=0; n<nCplx; ++n) {
    out[n_zeta - 1 - n] = out[n];
}
```

The output data is compared in the following figure with the equivalent result from using a `c2r` transform:

![R of magnetic axis: c2r DFT vs. REDFT01](img/axis_R.png)

The violet curve labelled 'DFT' is the output from the `c2r` DFT applied in the first half of this example.
The green curve labelled 'REDFT01' is the output from the `REDFT01` applied in the second half of this example.
Note that the locations of the `REDFT01` output needed to be shifted to the right by half an index to yield overlapping curves.

The full example can be found in [`src/app_magn_axis_stellsym.c`](src/app_magn_axis_stellsym.c).

### Geometry of a Flux Surface in a Stellarator
The real-space geometry of the flux surface is a general toroidal surface
which is conveniently parameterized using poloidal and toroidal angle-like coordinates.
The surface geometry is then given via a two-dimensional DFT:

![full Fourier series for flux surface](eqn/flux_surface.png)

where (θ,ζ) are the poloidal and toroidal angle-like variables in radians, respectively.
For a five-fold symmetric stellarator like W7-X, ζ ranges from 0 to 5 * 2π around the whole machine.
Conversely, the first (unique) toroidal segment of the geometry is contained within a range of ζ from 0 to 2π.
The poloidal angle-like coordinate θ ranges from 0 to 2π once the short way around the torus (wristband-like).

The two-dimensional DFT above can be simplified (shown here for *R* only)
by introducing rectangular two-dimensional Fourier coefficient matrices:

![full 2d FC matrix for flux surface](eqn/vmec_2d_fc.png)

In this example the geometry of the flux surface is evaluated
on a regular grid in θ and ζ:

![grid in zeta](eqn/grid_zeta.png)

![grid in theta](eqn/grid_theta.png)

where *l* ranges from `0` to `n_ζ-1` and *k* ranges from `0` to `n_θ-1`.

Inserting this into the Fourier series for, e.g., the *R* coordinate, leads to a DFT:

![flux surface R via DFT](eqn/flux_surface_dft.png)

The *R* coordinate is a real-valued quantity
which implies that a two-dimensional `c2r` DFT provided by FFTW can be used to
perform the backward transform in order to evaluate the flux surface geometry.
One further issue consists in the fact that the definition
of the Fourier geometry employed in VMEC uses an angle argument *(m θ - n ζ)*
where the two-dimensional DFT written out above uses *(m θ + n ζ)*.
The Fourier coefficients coming from VMEC
have to be inserted into the positions corresponding to *(-n)*
in the input array given to FFTW in order to resolve this pecularity. 

The next part concers the code to copy over the Fourier coefficients
staring from the VMEC output in arrays `rmnc` and `zmns` into the input to FFTW.

The index in the FFTW input corresponding to given mode numbers `n` and `m`
can be computed as follows:
 
```C
if (n <= 0) {
    idx_in = -n * nyq_pol + m;
} else {
    idx_in = (n_zeta - n) * nyq_pol + m;
}
```

This code is abbreviated in the following as `/* compute idx_in */`.

For the first coefficients with `m=0`, 
the copying is performed as follows:

```C
int idx_vmec = 0;
int m = 0;
for (int n = 0; n <= ntor; ++n) {

    /* compute idx_in */

    in_R[idx_in] =     rmnc[idx_vmec];
    in_Z[idx_in] = I * zmns[idx_vmec];

    idx_vmec++;
}
```

Note that VMEC weights the coefficients for `m=0` with positive `n ζ`
in contrast to the following coefficients with `m>0`.
The cosine-parity quantities like *R* are not influence by this,
but the sine-parity quantities like *Z* need to get the sign of their value flipped
due to the odd symmetry `sin(-n zeta) = -sin(n zeta)`.

The second part of the copying adresses the Fourier coefficients
with `m>0`:

```C
for (m = 1; m < mpol; ++m) {
    for (int n = -ntor; n <= ntor; ++n) {

        /* compute idx_in */

        in_R[idx_in] =  0.5 *     rmnc[idx_vmec];
        in_Z[idx_in] = -0.5 * I * zmns[idx_vmec];

        idx_vmec++;
    }
}
```

Here, the coefficients are scaled by a factor of `0.5`.
FFTW implies that coefficientss for `m<0` were present in the logically equivalent DFT input,
which is not the case for VMEC output. Thus, they are (in this case errornously)
scaled internally by a factor of `2`, which needs to be counteracted by a scaling factor of `0.5` in the input.

Assuming [stellarator symmetry](https://doi.org/10.1016/S0167-2789(97)00216-9),
half of the Fourier coefficients can be omitted and the transform reduces to the two-dimensional IDCT and IDST:

![stellarator-symmetric Fourier series for flux surface](eqn/flux_surface_stellsym.png)

This is currently resolved by simply omitting the stellarator-asymmetric terms
in the input to FFTW.

The full example code can be found in [`src/app_flux_surface.c`](src/app_flux_surface.c).
A Python script is provided to plot the resulting real-space geometry in [`src/plot_lcfs_realspace.py`](src/plot_lcfs_realspace.py).

## Allocation of arrays
Throughout this example collection, the proposed convenience wrapper functions provided by FFTW for allocating real- and complex-valued arrays are used:
```C
int n = 32;
int nOut = n/2+1;
double *in = fftw_alloc_real(n);
fftw_complex *out = fftw_alloc_complex(nOut);
```
where `n` is the real-space size of the DFT and `nOut` is the number of Fourier coefficients resulting from a `r2c` DFT.
The corresponding "raw" allocation code would look like this:
```C
double *in = (double*) fftw_malloc(sizeof(double) * n);
fftw_complex *out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nOut);
```
Note that above code is equivalent to the standard C way of allocating memory using `malloc`:
```C
double *in = (double*) malloc(sizeof(double) * n);
fftw_complex *out = (fftw_complex*) malloc(sizeof(fftw_complex) * nOut);
```
except that the FFTW routines ensure proper memory alignment for exploiting SIMD instructions of modern CPUs.

## Utility functions
In order to keep the examples short, a separate header file [`src/util.h`](src/util.h) is provided.
It contains methods to operate on one- and two-dimensional arrays (the latter stored in row-major order)
of real (`double`) and complex (`fftw_complex`) numbers.
The following operations are supported:
 * fill with random numbers between 0 and 1: e.g. `fill_random_1d_cplx`
 * element-wise check for approximate equality: e.g. `compare_1d_real`
 * write into a text file: e.g. `dump_1d_real`
 * compute generalized DFT (shift and non-linear phase): `dft_1d_cplx`

