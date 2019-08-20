// alloc.cc
/* Routines for allocation and deallocation of dynamically allocated memory. */
#include <complex>
#include <iostream>
#include "alloc.h" // Include header for consistency check
using std::complex;


// Allocation of 3D arrays ********************

double*** Alloc3D_d(const int Dim0Len, const int Dim1Len, const int Dim2Len)
{
    /* Allocates memory contiguously for a 3d array (size Dim0Len*Dim1Len*Dim2Len) of 
       doubles, in row-major layout (i.e. the last index changes fastest). */
    
    // The array is pointed to by a double***.
    double*** Array = new double** [Dim0Len];
    
    // Each row gets an array of pointers.
    for (int i=0; i<Dim0Len; ++i) Array[i] = new double* [Dim1Len];
    
    // Memory for the entire array is reserved at the first element.
    Array[0][0] = new double [Dim0Len*Dim1Len*Dim2Len];
    
    // The pointers for rows and columns are assigned to the appropriate location in the grid.
    // First row:
    for (int j=1; j<Dim1Len; ++j)
    {
        Array[0][j] = &Array[0][0][j*Dim2Len];
    }
    
    // Other rows:
    for (int i=1; i<Dim0Len; ++i)
    {
        for (int j=0; j<Dim1Len; ++j)
        {
            Array[i][j] = &Array[0][0][i*Dim1Len*Dim2Len + j*Dim2Len];
        }
    }
    
    // The three-dimensional array is returned.
    return Array;
}
void Dealloc3D(double const *const *const *const Array, const int Dim0Len)
{
    /* Deallocates memory allocated by Allocate3dArray(), i.e. for a contiguously stored 
       3d array (size Dim0Len*Dim1Len*Dim2Len) of doubles in row-major layout (i.e. the 
       last index changes fastest). The memory is deallocated in the order opposite to 
       its allocation. */
    delete [] Array[0][0];
    for (int i=0; i<Dim0Len; ++i) delete [] Array[i];
    delete Array;
}

float*** Alloc3D_f(const int Dim0Len, const int Dim1Len, const int Dim2Len)
{
    /* Allocates memory contiguously for a 3d array (size Dim0Len*Dim1Len*Dim2Len) of 
       floats, in row-major layout (i.e. the last index changes fastest). */
    
    // The array is pointed to by a double***.
    float*** Array = new float** [Dim0Len];
    
    // Each row gets an array of pointers.
    for (int i=0; i<Dim0Len; ++i) Array[i] = new float* [Dim1Len];
    
    // Memory for the entire array is reserved at the first element.
    Array[0][0] = new float [Dim0Len*Dim1Len*Dim2Len];
    
    // The pointers for rows and columns are assigned to the appropriate location in the grid.
    // First row:
    for (int j=1; j<Dim1Len; ++j)
    {
        Array[0][j] = &Array[0][0][j*Dim2Len];
    }
    
    // Other rows:
    for (int i=1; i<Dim0Len; ++i)
    {
        for (int j=0; j<Dim1Len; ++j)
        {
            Array[i][j] = &Array[0][0][i*Dim1Len*Dim2Len + j*Dim2Len];
        }
    }
    
    // The three-dimensional array is returned.
    return Array;
}
void Dealloc3D(float const *const *const *const Array, const int Dim0Len)
{
    /* Deallocates memory allocated by Alloc3d_f(), i.e. for a contiguously stored 
       3d array (size Dim0Len*Dim1Len*Dim2Len) of floats in row-major layout (i.e. the 
       last index changes fastest). The memory is deallocated in the order opposite to 
       its allocation. */
    delete [] Array[0][0];
    for (int i=0; i<Dim0Len; ++i) delete [] Array[i];
    delete Array;
}

complex<double>*** Alloc3D_z(const int Dim0Len, const int Dim1Len, const int Dim2Len)
{
    /* Allocates memory contiguously for a 3d array (size Dim0Len*Dim1Len*Dim2Len) of 
       doubles, in row-major layout (i.e. the last index changes fastest). */
    
    // The array is pointed to by a double***.
    complex<double>*** Array = new complex<double>** [Dim0Len];
    
    // Each row gets an array of pointers.
    for (int i=0; i<Dim0Len; ++i) Array[i] = new complex<double>* [Dim1Len];
    
    // Memory for the entire array is reserved at the first element.
    Array[0][0] = new complex<double> [Dim0Len*Dim1Len*Dim2Len];
    
    // The pointers for rows and columns are assigned to the appropriate location in the grid.
    // First row:
    for (int j=1; j<Dim1Len; ++j)
    {
        Array[0][j] = &Array[0][0][j*Dim2Len];
    }
    
    // Other rows:
    for (int i=1; i<Dim0Len; ++i)
    {
        for (int j=0; j<Dim1Len; ++j)
        {
            Array[i][j] = &Array[0][0][i*Dim1Len*Dim2Len + j*Dim2Len];
        }
    }
    
    // The three-dimensional array is returned.
    return Array;
}
void Dealloc3D(complex<double> const *const *const *const Array, const int Dim0Len)
{
    /* Deallocates memory allocated by Allocate3dArray(), i.e. for a contiguously stored 
       3d array (size Dim0Len*Dim1Len*Dim2Len) of doubles in row-major layout (i.e. the 
       last index changes fastest). The memory is deallocated in the order opposite to 
       its allocation. */
    delete [] Array[0][0];
    for (int i=0; i<Dim0Len; ++i) delete [] Array[i];
    delete Array;
}

complex<float>*** Alloc3D_c(const int Dim0Len, const int Dim1Len, const int Dim2Len)
{
    /* Allocates memory contiguously for a 3d array (size Dim0Len*Dim1Len*Dim2Len) of 
       doubles, in row-major layout (i.e. the last index changes fastest). */
    
    // The array is pointed to by a double***.
    complex<float>*** Array = new complex<float>** [Dim0Len];
    
    // Each row gets an array of pointers.
    for (int i=0; i<Dim0Len; ++i) Array[i] = new complex<float>* [Dim1Len];
    
    // Memory for the entire array is reserved at the first element.
    Array[0][0] = new complex<float> [Dim0Len*Dim1Len*Dim2Len];
    
    // The pointers for rows and columns are assigned to the appropriate location in the grid.
    // First row:
    for (int j=1; j<Dim1Len; ++j)
    {
        Array[0][j] = &Array[0][0][j*Dim2Len];
    }
    
    // Other rows:
    for (int i=1; i<Dim0Len; ++i)
    {
        for (int j=0; j<Dim1Len; ++j)
        {
            Array[i][j] = &Array[0][0][i*Dim1Len*Dim2Len + j*Dim2Len];
        }
    }
    
    // The three-dimensional array is returned.
    return Array;
}
void Dealloc3D(complex<float> const *const *const *const Array, const int Dim0Len)
{
    /* Deallocates memory allocated by Allocate3dArray(), i.e. for a contiguously stored 
       3d array (size Dim0Len*Dim1Len*Dim2Len) of doubles in row-major layout (i.e. the 
       last index changes fastest). The memory is deallocated in the order opposite to 
       its allocation. */
    delete [] Array[0][0];
    for (int i=0; i<Dim0Len; ++i) delete [] Array[i];
    delete Array;
}


// Allocation of 2D arrays ********************

double** Alloc2D_d(const int NumRows, const int NumCols)
{
    /* Allocates memory contiguously for a 2d array of size NumRows*NumCols. */
    double** Array = new double* [NumRows];
    
    // Memory for the entire array is reserved at the first element.
    Array[0] = new double [NumRows*NumCols];
    // The pointers for the subsequent elements of Array are assigned to the appropriate location.
    for (int i=1; i<NumRows; ++i) Array[i] = &Array[0][i*NumCols];
    
    // Two-dimensional array is returned.
    return Array;
}
void Dealloc2D(double const *const *const Array)
{
    /* Deallocates memory for a contiguous 1d grid of arrays (i.e. a 2d array) of type 
       "double". Memory is released in the order opposite to its allocation. */
    delete [] Array[0];
    delete [] Array;
}

float** Alloc2D_f(const int NumRows, const int NumCols)
{
    /* Allocates memory contiguously for a 2d array of size NumRows*NumCols. */
    float** Array = new float* [NumRows];
    
    // Memory for the entire array is reserved at the first element.
    Array[0] = new float [NumRows*NumCols];
    // The pointers for the subsequent elements of Array are assigned to the appropriate location.
    for (int i=1; i<NumRows; ++i) Array[i] = &Array[0][i*NumCols];
    
    // Two-dimensional array is returned.
    return Array;
}
void Dealloc2D(float const *const *const Array)
{
    /* Deallocates memory for a contiguous 2d array of type "float". Memory is released 
    in the order opposite to its allocation. */
    delete [] Array[0];
    delete [] Array;
}

std::complex<double>** Alloc2D_z(const int NumRows, const int NumCols)
{
    /* Allocates memory contiguously in row-major layout for a 2d array of size NumRows*NumCols. */
    std::complex<double>** Array = new std::complex<double>* [NumRows];
    
    // Memory for the entire array is reserved at the first element.
    Array[0] = new std::complex<double> [NumRows*NumCols];
    // The pointers for the subsequent elements of Array are assigned to the appropriate location.
    for (int i=1; i<NumRows; ++i) Array[i] = &Array[0][i*NumCols];
    
    // Two-dimensional array is returned.
    return Array;
}
void Dealloc2D(std::complex<double> const *const *const Array)
{
    /* Deallocates memory for a contiguous 1d grid of arrays (i.e. a 2d array) of type 
       "double". Memory is released in the order opposite to its allocation. */
    delete [] Array[0];
    delete [] Array;
}

std::complex<float>** Alloc2D_c(const int NumRows, const int NumCols)
{
    /* Allocates memory contiguously in row-major layout for a 2d array of size NumRows*NumCols. */
    std::complex<float>** Array = new std::complex<float>* [NumRows];
    
    // Memory for the entire array is reserved at the first element.
    Array[0] = new std::complex<float> [NumRows*NumCols];
    // The pointers for the subsequent elements of Array are assigned to the appropriate location.
    for (int i=1; i<NumRows; ++i) Array[i] = &Array[0][i*NumCols];
    
    // Two-dimensional array is returned.
    return Array;
}
void Dealloc2D(std::complex<float> const *const *const Array)
{
    /* Deallocates memory for a contiguous 1d grid of arrays (i.e. a 2d array) of type 
       "double". Memory is released in the order opposite to its allocation. */
    delete [] Array[0];
    delete [] Array;
}
