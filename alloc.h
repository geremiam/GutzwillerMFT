// alloc.h
/* Module containing routines for allocating and deallocating arrays of varying 
   dimensionality. */
#ifndef ALLOC_H
#define ALLOC_H

#include <complex> // Requires the "complex" library.


// Allocation of 3D arrays ********************

double*** Alloc3D_d(const int Dim0Len, const int Dim1Len, const int Dim2Len);
void Dealloc3D(double const *const *const *const Array, const int Dim0Len);

float***  Alloc3D_f(const int Dim0Len, const int Dim1Len, const int Dim2Len);
void Dealloc3D(float const *const *const *const Array, const int Dim0Len);

std::complex<double>*** Alloc3D_z(const int Dim0Len, const int Dim1Len, const int Dim2Len);
void Dealloc3D(std::complex<double> const *const *const *const Array, const int Dim0Len);

std::complex<float>*** Alloc3D_c(const int Dim0Len, const int Dim1Len, const int Dim2Len);
void Dealloc3D(std::complex<float> const *const *const *const Array, const int Dim0Len);


// Allocation of 2D arrays ********************

double** Alloc2D_d(const int NumRows, const int NumCols);
void Dealloc2D(double const *const *const Array);

float** Alloc2D_f(const int NumRows, const int NumCols);
void Dealloc2D(float const *const *const Array);

std::complex<double>** Alloc2D_z(const int NumRows, const int NumCols);
void Dealloc2D(std::complex<double> const *const *const Array);

std::complex<float>** Alloc2D_c(const int NumRows, const int NumCols);
void Dealloc2D(std::complex<float> const *const *const Array);


#endif
