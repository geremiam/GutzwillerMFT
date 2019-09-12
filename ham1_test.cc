// ham1_test.cc
#include <iostream>
#include <iomanip> // For the function std::setprecision()
#include <complex> // For complex numbers
#include <cmath> // For the constant M_PI
#include "IO.h"
#include "ham1.h"
using std::complex;

void MFs_t_test_1()
{
    // Construct without arguments
    MFs_t MFs;
    
    MFs.chi_s = {1., 2.};
    
    std::cout << MFs << std::endl;
}
void MFs_t_test_2()
{
    // Construct with arguments
    MFs_t MFs(1., 2., 3., 4.);
    
    std::cout << MFs << std::endl;
}
void MFs_t_test_3()
{
    // Construct an array
    
    complex<double>*const array = new complex<double> [4];
    
    array[0] = {0.,1.};
    array[1] = {0.,2.};
    array[2] = {0.,3.};
    array[3] = {0.,4.};
    
    MFs_t MFs(array);
    
    std::cout << MFs << std::endl;
    
    // Test the method MFs_to_array 
    complex<double> x [4] = {0.};
    MFs.MFs_to_array(x);
    PrintVector(4, x, std::cout);
    
    delete [] array;
}
void MFs_t_test_4()
{ // Testing assignment and copying
    // Construct with arguments
    MFs_t MFs(-1., 2., -3., 4.);
    std::cout << "MFs:\t" << MFs << std::endl;
    
    MFs_t MF1(MFs);
    std::cout << "MF1:\t" << MF1 << std::endl;
    
    MFs_t MF2 = MFs;
    std::cout << "MF2:\t" << MF2 << std::endl;
}
void MFs_t_test_5()
{ // Testing assignment and copying
    // Construct with arguments
    MFs_t MFs(1., 2., 3., 4.);
    std::cout <<  "MFs = \t" <<  MFs << std::endl;
    std::cout << "-MFs = \t" << -MFs << std::endl;
    std::cout <<  "MFs = \t" <<  MFs << std::endl;
    std::cout << "MFs*0.5 = \t" << MFs*0.5 << std::endl;
    std::cout <<  "MFs = \t" <<  MFs << std::endl;
    
    MFs_t MF2 = MFs+MFs;
    std::cout << "MF2 = \t" << MF2 << std::endl;
    std::cout <<  "MFs = \t" <<  MFs << std::endl;
    
    MFs_t MF3 = (MFs*1.2 + MF2)*0.1;
    std::cout <<  "MFs = \t" <<  MF3 << std::endl;
    const bool is_small = MF3.check_bound(1.3);
    std::cout << is_small << std::endl;
}






int main()
{
    MFs_t_test_5();
    return 0;
}
