// IO.cc
/* Input and output in different formats */
#include <complex>
#include <ostream>
#include "IO.h"

void PrintMatrix(const int num_rows, const int num_cols, const int*const mat, std::ostream& out)
{
    /* Prints a matrix held as a 1D array in row-major storage. */
    for (int index=0; index<num_rows*num_cols; ++index)
    {
        int row = index/num_cols;
        int col = index%num_cols;
        
        out << mat[index];
        if (col==num_cols-1)
            out << '\n';
        else
            out << " ";
    }
}

void PrintMatrix(const int num_rows, const int num_cols, const int*const*const mat, std::ostream& out)
{
    /* Prints a matrix held as a 2D array in row-major storage. */
    for (int row=0; row<num_rows; ++row)
    {
        for (int col=0; col<num_cols; ++col)
        {
            out << mat[row][col];
            if (col!=num_cols-1) out << " ";
        }
        out << "\n";
    }
}

void PrintMatrix(const int num_rows, const int num_cols, const double*const mat, std::ostream& out)
{
    /* Prints a matrix held as a 1D array in row-major storage. */
    for (int index=0; index<num_rows*num_cols; ++index)
    {
        int row = index/num_cols;
        int col = index%num_cols;
        
        out << mat[index];
        if (col==num_cols-1)
            out << '\n';
        else
            out << " ";
    }
}

void PrintMatrix(const int num_rows, const int num_cols, const double*const*const mat, std::ostream& out)
{
    /* Prints a matrix held as a 2D array in row-major storage. */
    for (int row=0; row<num_rows; ++row)
    {
        for (int col=0; col<num_cols; ++col)
        {
            out << mat[row][col];
            if (col!=num_cols-1) out << " ";
        }
        out << "\n";
    }
}

void PrintMatrix(const int num_rows, const int num_cols, const std::complex<double>*const mat, std::ostream& out)
{
    /* Prints a matrix held as a 1D array in row-major storage. */
    for (int index=0; index<num_rows*num_cols; ++index)
    {
        int row = index/num_cols;
        int col = index%num_cols;
        
        out << mat[index];
        if (col==num_cols-1)
            out << '\n';
        else
            out << "\t";
    }
}

void PrintMatrix(const int num_rows, const int num_cols, const std::complex<double>*const*const mat, std::ostream& out)
{
    /* Prints a matrix held as a 2D array in row-major storage. */
    for (int row=0; row<num_rows; ++row)
    {
        for (int col=0; col<num_cols; ++col)
        {
            out << mat[row][col];
            if (col!=num_cols-1) out << "\t";
        }
        out << "\n";
    }
}


void PrintVector(const int len, const double*const vec, std::ostream& out)
{
    /* Prints a vector held in a 1D array. */
    out << "{";
    for (int index=0; index<len; ++index)
    {
        out << vec[index];
        if (index!=len-1) out << " ";
    }
    out << "}" << std::endl;
}

void PrintVector(const int len, const float*const vec, std::ostream& out)
{
    /* Prints a vector held in a 1D array. */
    out << "{";
    for (int index=0; index<len; ++index)
    {
        out << vec[index];
        if (index!=len-1) out << " ";
    }
    out << "}" << std::endl;
}

void PrintVector(const int len, const std::complex<double>*const vec, std::ostream& out)
{
    /* Prints a vector held in a 1D array. */
    out << "{";
    for (int index=0; index<len; ++index)
    {
        out << vec[index];
        if (index!=len-1) out << "\t";
    }
    out << "}" << std::endl;
}

void PrintVector(const int len, const std::complex<float>*const vec, std::ostream& out)
{
    /* Prints a vector held in a 1D array. */
    out << "{";
    for (int index=0; index<len; ++index)
    {
        out << vec[index];
        if (index!=len-1) out << "\t";
    }
    out << "}" << std::endl;
}
