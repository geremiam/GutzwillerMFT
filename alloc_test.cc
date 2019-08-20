// alloc_test.cc
// Testing suite for the module alloc_dealloc.
#include <iostream>
#include <complex>
#include "alloc.h"

void Alloc3D_d_TEST()
{
    // Unit test for the routine Alloc3d_d(): we construct a 3*5*2 array.
    
    // Allocate the array
    const int Dim0Len = 3;
    const int Dim1Len = 5;
    const int Dim2Len = 2;
    double*** A = Alloc3D_d(Dim0Len, Dim1Len, Dim2Len);
    
    // Index values are assigned to the array, starting from the very first entry
    for (int i=0; i<Dim0Len*Dim1Len*Dim2Len; ++i)
    {
        A[0][0][i] = (double)(i);
    }
    
    // We check that the indices have been assigned as expected:
    for (int i=0; i<Dim0Len; ++i)
    {
        for (int j=0; j<Dim1Len; ++j)
        {
            for (int k=0; k<Dim2Len; ++k)
                std::cout << A[i][j][k] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
    
    // Memory deallocation
    Dealloc3D(A, Dim0Len);
}

void Alloc3D_f_TEST()
{
    // Unit test for the routine Alloc3d_f(): we construct a 3*5*2 array.
    
    // Allocate the array
    const int Dim0Len = 3;
    const int Dim1Len = 5;
    const int Dim2Len = 2;
    float***  A = Alloc3D_f(Dim0Len, Dim1Len, Dim2Len);
    
    // Index values are assigned to the array, starting from the very first entry
    for (int i=0; i<Dim0Len*Dim1Len*Dim2Len; ++i)
    {
        A[0][0][i] = (float)(i);
    }
    
    // We check that the indices have been assigned as expected:
    for (int i=0; i<Dim0Len; ++i)
    {
        for (int j=0; j<Dim1Len; ++j)
        {
            for (int k=0; k<Dim2Len; ++k)
                std::cout << A[i][j][k] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
    
    // Memory deallocation
    Dealloc3D(A, Dim0Len);
}

void Alloc3D_z_TEST()
{
    // Unit test for the routine Alloc3d_d(): we construct a 3*5*2 array.
    
    // Allocate the array
    const int Dim0Len = 3;
    const int Dim1Len = 5;
    const int Dim2Len = 2;
    std::complex<double>*** A = Alloc3D_z(Dim0Len, Dim1Len, Dim2Len);
    
    // Index values are assigned to the array, starting from the very first entry
    for (int i=0; i<Dim0Len*Dim1Len*Dim2Len; ++i)
    {
        A[0][0][i] = {(double)(i), 0.1*(double)(i)};
    }
    
    // We check that the indices have been assigned as expected:
    for (int i=0; i<Dim0Len; ++i)
    {
        for (int j=0; j<Dim1Len; ++j)
        {
            for (int k=0; k<Dim2Len; ++k)
                std::cout << A[i][j][k] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
    
    // Memory deallocation
    Dealloc3D(A, Dim0Len);
}

void Alloc3D_c_TEST()
{
    // Unit test for the routine Alloc3d_d(): we construct a 3*5*2 array.
    
    // Allocate the array
    const int Dim0Len = 3;
    const int Dim1Len = 5;
    const int Dim2Len = 2;
    std::complex<float>*** A = Alloc3D_c(Dim0Len, Dim1Len, Dim2Len);
    
    // Index values are assigned to the array, starting from the very first entry
    for (int i=0; i<Dim0Len*Dim1Len*Dim2Len; ++i)
    {
        A[0][0][i] = {(float)(i), 0.1f*(float)(i)};
    }
    
    // We check that the indices have been assigned as expected:
    for (int i=0; i<Dim0Len; ++i)
    {
        for (int j=0; j<Dim1Len; ++j)
        {
            for (int k=0; k<Dim2Len; ++k)
                std::cout << A[i][j][k] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl << std::endl;
    }
    
    // Memory deallocation
    Dealloc3D(A, Dim0Len);
}


void Alloc2D_d_TEST()
{
    /* Unit test for the routine Alloc2d_d: the indexing is checked for a 2*2
       case. */
    const int NumRows = 2;
    const int NumCols = 2;
    double** A = Alloc2D_d(NumRows, NumCols);
    
    // Values 0 to 3 are assigned to the array associated to the first element.
    for (int i=0; i<NumRows*NumCols; ++i) A[0][i] = (double)(i);
    
    // Indexing check:
    for (int i=0; i<NumRows; ++i)
    {
        for (int j=0; j<NumCols; ++j)
            std::cout << A[i][j] << " ";
        std::cout << std::endl;
    }
    
    // Memory deallocation
    Dealloc2D(A);
}

void Alloc2D_f_TEST()
{
    /* Unit test for the routine Alloc2d_f: the indexing is checked for a 2*2
       case. */
    
    const int NumRows = 2;
    const int NumCols = 8;
    float**  A = Alloc2D_f(NumRows, NumCols);
    // Values are assigned to the array associated to the first element.
    for (int i=0; i<NumRows*NumCols; ++i) A[0][i] = (double)(i);
    
    // Indexing check:
    for (int i=0; i<NumRows; ++i)
    {
        for (int j=0; j<NumCols; ++j)
            std::cout << A[i][j] << " ";
        std::cout << std::endl;
    }
    
    // Memory deallocation
    Dealloc2D(A);
}

void Alloc2D_z_TEST()
{
    /* Unit test for the routine Alloc2d_z: the indexing is checked for a 2*2
       case. */
    
    const int NumRows = 2;
    const int NumCols = 8;
    std::complex<double>**  A = Alloc2D_z(NumRows, NumCols);
    // Values are assigned to the array associated to the first element.
    for (int i=0; i<NumRows*NumCols; ++i) A[0][i] = {(double)(i), -0.5*(double)(i)};
    
    // Indexing check:
    for (int i=0; i<NumRows; ++i)
    {
        for (int j=0; j<NumCols; ++j)
            std::cout << A[i][j] << " ";
        std::cout << std::endl;
    }
    
    // Memory deallocation
    Dealloc2D(A);
}

void Alloc2D_c_TEST()
{
    /* Unit test for the routine Alloc2d_z: the indexing is checked for a 2*2
       case. */
    
    const int NumRows = 2;
    const int NumCols = 8;
    std::complex<float>**  A = Alloc2D_c(NumRows, NumCols);
    // Values are assigned to the array associated to the first element.
    for (int i=0; i<NumRows*NumCols; ++i) A[0][i] = {(float)(i), -0.5f*(float)(i)};
    
    // Indexing check:
    for (int i=0; i<NumRows; ++i)
    {
        for (int j=0; j<NumCols; ++j)
            std::cout << A[i][j] << " ";
        std::cout << std::endl;
    }
    
    // Memory deallocation
    Dealloc2D(A);
}


int main()
{
    Alloc3D_d_TEST(); std::cout << "_-*-_-*-_-*-_-*-_-*-_" << std::endl;
    
    Alloc3D_f_TEST(); std::cout << "_-*-_-*-_-*-_-*-_-*-_" << std::endl;
    
    Alloc3D_z_TEST(); std::cout << "_-*-_-*-_-*-_-*-_-*-_" << std::endl;
    
    Alloc3D_c_TEST(); std::cout << "_-*-_-*-_-*-_-*-_-*-_" << std::endl;
    
    Alloc2D_d_TEST(); std::cout << "_-*-_-*-_-*-_-*-_-*-_" << std::endl;
    
    Alloc2D_f_TEST(); std::cout << "_-*-_-*-_-*-_-*-_-*-_" << std::endl;
    
    Alloc2D_z_TEST(); std::cout << "_-*-_-*-_-*-_-*-_-*-_" << std::endl;
    
    Alloc2D_c_TEST(); std::cout << "_-*-_-*-_-*-_-*-_-*-_" << std::endl;
    
    return 0;
}
