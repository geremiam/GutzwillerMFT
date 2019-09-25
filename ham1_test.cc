// ham1_test.cc
#include <iostream>
#include <iomanip> // For the function std::setprecision()
#include <complex> // For complex numbers
#include <cmath> // For the constant M_PI
//#include "IO.h"
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
    std::cout <<  "MFs = \t" <<  MFs << std::endl << std::endl;
    
    std::cout << "-MFs = \t" << -MFs << std::endl;
    std::cout <<  "MFs = \t" <<  MFs << std::endl << std::endl;
    
    std::cout << "MFs*0.5 = \t" << MFs*0.5 << std::endl;
    std::cout <<  "MFs = \t" <<  MFs << std::endl << std::endl;
    
    MFs_t MF2 = MFs+MFs;
    std::cout << "MF2 = \t" << MF2 << std::endl;
    std::cout <<  "MFs = \t" <<  MFs << std::endl << std::endl;
    
    MFs_t MF3 = MFs - MFs*2.;
    std::cout << "MF3 = \t" << MF3 << std::endl;
    std::cout <<  "MFs = \t" <<  MFs << std::endl << std::endl;
    
    MFs_t MF4 = (MFs*1.2 - MF2)*0.1;
    std::cout <<  "MF4 = \t" <<  MF4 << std::endl;
    
    //const bool is_small = MF3.check_bound(1.3);
    //std::cout << is_small << std::endl;
}

// **************************************************************************************
void FPparams_t_test_1()
{
    const double tol = 1.e-6;
    const int loops_lim = 300;
    const int mixing_vals_len = 4;
    const int counter_vals [mixing_vals_len] = {10, 20, 30, 40};
    const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6};
    FPparams_t FPparams(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    
    for (int i=0; i<50; ++i)
        std::cout << i << "\t" << FPparams.mixratio(i) << std::endl;
}
void FPparams_t_test_2()
{
    // Testing copy constructor
    const double tol = 1.e-6;
    const int loops_lim = 300;
    const int mixing_vals_len = 4;
    const int counter_vals [mixing_vals_len] = {10, 20, 30, 40};
    const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6};
    
    FPparams_t FPparam1(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    std::cout << "FPparam1:\n" << FPparam1 << std::endl << std::endl;
    
    
    FPparams_t FPparam2(FPparam1);
    std::cout << "FPparam2:\n" << FPparam2 << std::endl << std::endl;
}


// **************************************************************************************
void ham1_t_test_1()
{
    // Test for the diag method
    
    // We need to declare parameters for the fixed point algorithm.
    const double tol = 1.e-6;
    const int loops_lim = 300;
    const int mixing_vals_len = 4;
    const int counter_vals [mixing_vals_len] = {10, 20, 30, 40};
    const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6};
    FPparams_t FPparams(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    
    // Declare an object to hold the initial values for the MFs.
    MFs_t MFs_initial(0., 0., 0., 0.);
    
    // Declare an instance of the Hamiltonian
    ham1_t ham1(FPparams, MFs_initial, 500, 500);
    
    std::cout << ham1.MFs_ << std::endl;
    
    const double kx = 0.*M_PI;
    const double ky = 0.*M_PI;
    const double mu = 0.;
    
    double E_cal=0.;
    double u=0.;
    complex<double> v=0.;
    
    ham1.diag(kx, ky, mu, E_cal, u, v);
    
    std::cout << "E_cal = " << E_cal << ", u = " << u << ", v = " << v << std::endl;
}
void ham1_t_test_2()
{
    // Test of the chempot_utility method
    
    // We need to declare parameters for the fixed point algorithm.
    const double tol = 1.e-6;
    const int loops_lim = 300;
    const int mixing_vals_len = 4;
    const int counter_vals [mixing_vals_len] = {10, 20, 30, 40};
    const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6};
    FPparams_t FPparams(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    
    // Declare an object to hold the initial values for the MFs.
    MFs_t MFs_initial(0., 0., 0., 0.);
    
    // Declare an instance of the Hamiltonian
    ham1_t ham1(FPparams, MFs_initial, 500, 500);
    
    std::cout << ham1.MFs_ << std::endl;
    
    std::cout << ham1.chempot_utility(-1.) << std::endl;
}
void ham1_t_test_3()
{
    // Testfor the chempot method
    
    // We need to declare parameters for the fixed point algorithm.
    const double tol = 1.e-6;
    const int loops_lim = 300;
    const int mixing_vals_len = 4;
    const int counter_vals [mixing_vals_len] = {10, 20, 30, 40};
    const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6};
    FPparams_t FPparams(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    
    // Declare an object to hold the initial values for the MFs.
    MFs_t MFs_initial(0.1, 0.2, {0.,0.}, {0.,0.});
    
    // Declare an instance of the Hamiltonian
    ham1_t ham1(FPparams, MFs_initial, 500, 500);
    
    std::cout << ham1.MFs_ << std::endl;
    
    std::cout << "Calculating chemical potential.\n" << std::endl;
    std::cout << std::setprecision(15) << ham1.chempot(true) << std::endl;
}
void ham1_t_test_4()
{
    // Test for the compute_MFs() method
    
    // We need to declare parameters for the fixed point algorithm.
    const double tol = 1.e-6;
    const int loops_lim = 300;
    const int mixing_vals_len = 4;
    const int counter_vals [mixing_vals_len] = {10, 20, 30, 40};
    const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6};
    FPparams_t FPparams(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    
    // Declare an object to hold the initial values for the MFs.
    MFs_t MFs_initial(0.1, 0.2, 0., 0.);
    
    // Declare an instance of the Hamiltonian
    ham1_t ham1(FPparams, MFs_initial, 500, 500);
    
    double mu_out=0.;
    MFs_t output = ham1.compute_MFs(mu_out);
    
    std::cout << "output MFs: " << output << std::endl;
    std::cout << "mu_out = " << mu_out << std::endl;
}
void ham1_t_test_5()
{
    // Test for the FixedPoint() method
    
    // We need to declare parameters for the fixed point algorithm.
    const double tol = 1.e-6;
    const int loops_lim = 300;
    const int mixing_vals_len = 4;
    const int counter_vals [mixing_vals_len] = {10, 20, 30, 40};
    const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6};
    FPparams_t FPparams(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    
    // Declare an object to hold the initial values for the MFs.
    MFs_t MFs_initial(0.1, 0.2, {0.1,0.2}, {0.1,0.2});
    
    // Declare an instance of the Hamiltonian
    ham1_t ham1(FPparams, MFs_initial, 500, 500);
    
    bool output = ham1.FixedPoint(true);
    
    std::cout << "output = " << output << std::endl;
}



int main()
{
    ham1_t_test_5();
    return 0;
}
