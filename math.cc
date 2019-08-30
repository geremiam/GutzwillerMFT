// math_routines.cc
/*
Module containing various math functions or algorithms used to be used in other 
applications.
*/

#include <iostream> // Input/output to command line
#include <cmath> // Numerous math functions
#include <complex> // Includes complex numbers
#include "math.h" // Include header file for consistency check


double log_1p_exp(const double x)
{
    // Calculates the function log(1+exp(x)).
    // The function log1p takes care of x << -1.
    // If x >> 1, then answer is approximately x; use the conditional operator for this.
    return (x>log_max_double) ? x : log1p(exp(x));
}


double nF0(double const energy)
{
    /* Zero-temperature Fermi function. Overloaded for floats. */
    double ans = 1.;
    if (energy > 0.)
        ans = 0.;
    
    return ans;    
}
float  nF0(float const energy)
{
    /* Zero-temperature Fermi function. Overloaded for doubles. */
    float ans = 1.f;
    if (energy > 0.f)
        ans = 0.f;
    
    return ans;    
}
double nF(const double T, const double energy)
{
    // Fermi function for doubles. 
    return 1./( 1.+std::exp(energy/T) );
}
float nF(const float T, const float energy)
{
    // Fermi function for doubles. Use doubles to avoid overflow.
    const double arg = (double)(energy)/(double)(T);
    const double ans = 1./( 1. + std::exp(arg) );
    return (float)(ans);
}
