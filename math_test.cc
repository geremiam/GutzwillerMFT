// math_routines_test.cc
/* Testing suite for the module math_routines. */

#include <iostream> // For terminal input and output
#include <iomanip> // For setprecision
#include <complex>
#include "math.h" // Module's header

// Internal routines for specific tests

double nF_1(const double T, const double energy)
{
    // Fermi function for doubles. Use long doubles to avoid overflow.
    const long double arg = (long double)(energy)/(long double)(T);
    const long double ans = 1.l/( 1.l + std::exp(arg) );
    return (double)(ans);
}

double nF_2(const double T, const double energy)
{
    // Fermi function for doubles.
    return 1./( 1. + std::exp(energy/T) );
}

double nF_3(const double T, const double energy)
{
    return (energy>0.) ? std::exp(-energy/T)/( 1.+std::exp(-energy/T) ) : 1./( 1.+std::exp(energy/T) );
}

double nF_4(const double T, const double energy)
{
    const double arg = energy/T;
    double ans=0.;
    
    if (arg<-30.)
        ans = 1.-std::exp(arg);
    else if (arg>200.)
        ans = std::exp(-arg);
    else
        ans = 1./( 1. + std::exp(energy/T) );
    return ans;
}

void compare_Fermi_functions()
{
    // Test to compare different implementations of the Fermi function.
    /* Outcome: using icpc, there seems to be no clear advantage in using more clever
    implementations. */
    std::cout << std::setprecision(20);
    const int n = 21;
    const double lo = -40.; //-40.//700.
    const double hi = -25.; //-35.//710.
          double energies[n] = {0.};
    const double T = 1.;
    
    for (int i=0; i<n; ++i)
    {
        energies[i] = lo + (double)(i)*(hi-lo)/(double)(n-1);
        std::cout << energies[i] << "\n"
                  << nF_1(T, energies[i]) << "\n"
                  << nF_2(T, energies[i]) << "\n"
                  << nF_3(T, energies[i]) << "\n"
                  << nF_4(T, energies[i]) << "\n"
                  //<< std::exp(-energies[i]/T) << "\n"
                  << 1.-std::exp(energies[i]/T)  << "\n"
                  << std::endl;
    }
}

double foo(const double x)
{
    return x*x - 1.5;
}

void bisec_test()
{
    const double a_in = -1.23;
    const double b_in = 1.22;
    const double output = bisec(a_in, b_in, foo);
    
    std::cout << "output = " << output << std::endl;
}

// **************************************************************************************

void test_nF0()
{
    /* Test for the routine nF0 */
    
    const int ArrLen = 4;
    const double E [ArrLen] = {-1., -0.1, 0.2, 0.7};
    const float  e [ArrLen] = {-1., -0.1, 0.2, 0.7};
    
    const std::complex<double> z = {2., 3.};
    
    for (int i=0; i<ArrLen; ++i)
        std::cout << "nF0(" << E[i] << ") = " << z*nF0(E[i]) << "\t"
                  << "nF0(" << e[i] << ") = " << nF0(e[i]) << std::endl;
}

void test_nF()
{
    /* Test for the routine nF */
    
    const int ArrLen = 4;
    const double E [ArrLen] = {-1., -0.1, 0.2, 0.7};
    const float  e [ArrLen] = {-1., -0.1, 0.2, 0.7};
    
    const double T = 10.;
    const float  t = 10.f;
    
    for (int i=0; i<ArrLen; ++i)
        std::cout << "nF(T, " << E[i] << ") = " << nF(T, E[i]) << "\t"
                  << "nF(t, " << e[i] << ") = " << nF(t, e[i]) << std::endl;
    
}



int main()
{
    //compare_Fermi_functions();
    
    //test_nF0();
    
    //test_nF();
    
    bisec_test();
    
    return 0;
}
