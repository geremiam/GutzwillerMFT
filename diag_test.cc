// diag_test.cc
// Testing suite for the module diag_routines.
#include <iostream> // Input/output to command line
#include <complex> // For complex numbers
#include "diag.h"

void PrintMatrix(const int order, std::complex<double> const *const a);
void PrintVector(const int len, double const *const a);

void test_simple_zheev()
{
    const int order = 4;
    std::complex<double> a [order*order]= {{1,2},   {0,0},    {0,0},   {0,0},
                                           {3,4},   {5,6},    {0,0},   {0,0},
                                           {7,8},   {9,10},   {11,12}, {0,0},
                                           {13,14}, {15, 16}, {17,18}, {19,20}};
    std::complex<double> z [order*order];
    
    std::cout << "a =\n";
    PrintMatrix(order, a);
    
    double w [order]={0.};
    
    const bool OutputEvecs = true;
    simple_zheev(order, a, w, OutputEvecs, z);
    
    std::cout << "\nw = ";
    PrintVector(order, w);
    
    std::cout << "\na =\n";
    PrintMatrix(order, a);
    
    std::cout << "\nz =\n";
    PrintMatrix(order, z);
}

void test_simple_zheevr()
{
    const int order = 4;
    std::complex<double> a [order*order]= {{1,2},   {0,0},    {0,0},   {0,0},
                                           {3,4},   {5,6},    {0,0},   {0,0},
                                           {7,8},   {9,10},   {11,12}, {0,0},
                                           {13,14}, {15, 16}, {17,18}, {19,20}};
    std::complex<double> z [order*order];
    
    std::cout << "a =\n";
    PrintMatrix(order, a);
    
    double w [order]={0.};
    
    const bool OutputEvecs = true;
    simple_zheevr(order, a, w, OutputEvecs, z);
    
    std::cout << "\nw = ";
    PrintVector(order, w);
    
    std::cout << "\na =\n";
    PrintMatrix(order, a);
    
    std::cout << "\nz =\n";
    PrintMatrix(order, z);
}

void SpeedTest()
{
    const int order = 8;
    std::complex<double> a [order*order]= {{1,2},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
                                           {3,4},{5,6},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
                                           {7,8},{9,0},{1,2},{0,0},{0,0},{0,0},{0,0},{0,0},
                                           {3,4},{5,6},{7,8},{9,0},{0,0},{0,0},{0,0},{0,0},
                                           {1,2},{3,4},{5,6},{7,8},{9,0},{0,0},{0,0},{0,0},
                                           {1,2},{3,4},{5,6},{7,8},{9,0},{1,2},{0,0},{0,0},
                                           {3,4},{5,6},{7,8},{9,0},{1,2},{3,4},{5,6},{0,0},
                                           {7,8},{9,0},{1,2},{3,4},{5,6},{7,8},{9,0},{1,2}};
    std::complex<double> z [order*order];
    double w [order]={0.};
    
    //std::cout << "a =\n";
    //PrintMatrix(order, a);
    
    const bool OutputEvecs = true;
    for (int i=0; i<500000; ++i)
        simple_zheev(order, a, w, OutputEvecs, z);
}

void EqualityTest()
{
    const int order = 8; // Define a matrix
    std::complex<double> a1 [order*order]={{3,4},{0,0},  {0,0},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           {1,2},{3,4},  {0,0},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           
                                           {7,8},{0,0},  {3,4},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           {0,0},{5,6},  {1,2},{3,4},  {0,0},{0,0},  {0,0},{0,0},
                                           
                                           {1,2},{0,0},  {0,0},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           {0,0},{1,2},  {0,0},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           
                                           {0,0},{0,0},  {1,2},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           {0,0},{0,0},  {0,0},{1,2},  {0,0},{0,0},  {0,0},{0,0}};
    
    std::complex<double> a2 [order*order]={{3,4},{0,0},  {0,0},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           {1,2},{3,4},  {0,0},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           
                                           {7,8},{0,0},  {3,4},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           {0,0},{5,6},  {1,2},{3,4},  {0,0},{0,0},  {0,0},{0,0},
                                           
                                           {1,2},{0,0},  {0,0},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           {0,0},{1,2},  {0,0},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           
                                           {0,0},{0,0},  {1,2},{0,0},  {0,0},{0,0},  {0,0},{0,0},
                                           {0,0},{0,0},  {0,0},{1,2},  {0,0},{0,0},  {0,0},{0,0}};
    double w1 [order]={0.}; // Declare first set of evecs and evals arrays
    std::complex<double> z1 [order*order];
    
    double w2 [order]={0.}; // Declare second set of evecs and evals arrays
    std::complex<double> z2 [order*order];
    
    const bool OutputEvecs = true; // Do diagonalization
    simple_zheev(order, a1, w1, OutputEvecs, z1);
    simple_zheevr(order, a2
    , w2, OutputEvecs, z2);
    
    double d_w [order]={0.}; // Declare arrays for differences
    std::complex<double> d_z [order*order];
    
    for (int index=0; index<order; ++index)
        d_w[index] = w2[index] - w1[index];
    
    for (int index=0; index<order*order; ++index)
    {
        if ( abs(z2[index]+z1[index]) < abs(z2[index]-z1[index]))
            d_z[index] = z2[index] + z1[index];
        else
            d_z[index] = z2[index] - z1[index];
    }
    
    std::cout << "\nd_w = ";
    PrintVector(order, d_w);
    
    std::cout << "\nd_z =\n";
    PrintMatrix(order, d_z);
}

int main()
{
    test_simple_zheev();
    //std::cout << std::endl;
    test_simple_zheevr();
    //SpeedTest();
    //EqualityTest();
    return 0;
}


void PrintMatrix(const int order, std::complex<double> const *const a)
{
    for (int i=0; i<order*order; ++i)
        {
            std::cout << a[i] << "\t";
            if (i%order==(order-1)) std::cout << std::endl;
        }
}
void PrintVector(const int len, double const *const w)
{
    for (int i=0; i<len; ++i)
        std::cout << w[i] << " ";
    std::cout << std::endl;
}
