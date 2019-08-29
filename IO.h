// IO.h
/* Input and output in different formats */
#ifndef IO_H
#define IO_H

#include <complex> // Requires the "complex" library.

void PrintMatrix(const int num_rows, const int num_cols, const int*const       mat, std::ostream& out);
void PrintMatrix(const int num_rows, const int num_cols, const int*const*const mat, std::ostream& out);
void PrintMatrix(const int num_rows, const int num_cols, const double*const       mat, std::ostream& out);
void PrintMatrix(const int num_rows, const int num_cols, const double*const*const mat, std::ostream& out);
void PrintMatrix(const int num_rows, const int num_cols, const std::complex<double>*const       mat, std::ostream& out);
void PrintMatrix(const int num_rows, const int num_cols, const std::complex<double>*const*const mat, std::ostream& out);

void PrintVector(const int len, const double*const vec, std::ostream& out);
void PrintVector(const int len, const float*const  vec, std::ostream& out);
void PrintVector(const int len, const std::complex<double>*const vec, std::ostream& out);
void PrintVector(const int len, const std::complex<float>*const  vec, std::ostream& out);

#endif
