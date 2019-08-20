// IO_test.cc
/* Testing suite for the module misc. */
#include <iostream>
#include <complex>
#include "alloc.h" //Allocation and deallocation of arrays (must include complex beforehand)
#include "IO.h" // Include header file for consistency check

void test_PrintMatrix()
{
    /* Unit test for the (overloaded) function PrintMatrix() */
    const int num_rows = 2;
    const int num_cols = 2;
    
    double*const*const A = Alloc2D_d(num_rows, num_cols);
    for (int i=0; i<num_rows*num_cols; ++i) A[0][i] = (double)(i);
    PrintMatrix(num_rows, num_rows, A, std::cout);
    Dealloc2D(A);
    
    std::cout << std::endl;
    
    double*const B = new double [num_rows*num_cols];
    for (int i=0; i<num_rows*num_cols; ++i) B[i] = (double)(i*i);
    PrintMatrix(num_rows, num_rows, B, std::cout);
    delete [] B;
    
    std::cout << std::endl;
    
    std::complex<double>*const*const C = Alloc2D_z(num_rows, num_cols);
    for (int i=0; i<num_rows*num_cols; ++i) C[0][i] = {0.28*(double)(i), -1.23*(double)(i)};
    PrintMatrix(num_rows, num_rows, C, std::cout);
    Dealloc2D(C);
    
    std::cout << std::endl;
    
    std::complex<double>*const D = new std::complex<double> [num_rows*num_cols];
    for (int i=0; i<num_rows*num_cols; ++i) D[i] = {0.28*(double)(i), -1.23*(double)(i)};
    PrintMatrix(num_rows, num_rows, D, std::cout);
    delete [] D;
    
}

void test_PrintVector()
{
    /* Unit test for the (overloaded) routine PrintVector() */
    
    // Declare a test array
    const int len = 5;
    float arra [len] = {3., 6., 10., 5., -8.};
    PrintVector(len, arra, std::cout);
    
    std::complex<double> vec [len] = {{1.,-5.}, {2.,-1.}, {3.,0.}, {4.,6.}, {5.,-7.}};
    PrintVector(len,vec,std::cout);
    
}

int main()
{
    test_PrintMatrix();
    test_PrintVector();
    
    return 0;
}
