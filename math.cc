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


double bisec(const double a_in, const double b_in, double foo(const double), const bool show_output)
{
    // The chemical potential is found using the bisection method. Uses 100 times machine 
    // epsilon as relative tolerance.
    
    // Make sure that a < b.
    if ( !(a_in<b_in) )
      std::cout << "WARNING: function \"bisec\" requires a_in<b_in." << std::endl;
    // Make sure there is a zero between the two inputs.
    if ( !(foo(a_in)*foo(b_in)<0.) )
      std::cout << "WARNING: the function does not change sign between a_in and b_in." << std::endl;
    const bool increasing = (foo(a_in) < foo(b_in)); // Determine if the function is increasing
    
    // Starting values for a and b. Choose them slightly outside the energy range to be safe.
    double a = a_in;
    double b = b_in;
    
    int counter = 0;
    bool converged = false;
    
    do // Update via Newton's method
    {
      counter++; // Increment 'counter' by one
      
      const double midpoint = (a+b)/2.; // Get the midpoint between a and b
      // Find the image of the function at the midpoint
      const double image = foo(midpoint);
      
      if (show_output)
          std::cout << "a = " << a << ", b = " << b << ", b-a = " << b-a << "\t"
                    << "midpoint = " << midpoint << ", image = " << image;
      
      if (image==0.)
      {
        a = midpoint; // If image is 0, the loop terminates and the value is returned.
        b = midpoint;
      }
        
      else if (image<0.)
      {
        if (increasing)
          a = midpoint; // If 'image' is negative, 'midpoint' is the new 'a'
        else
          b = midpoint;
      }
      else if (image>0.)
      {
        if (increasing)
          b = midpoint; // If 'image' is positive, 'midpoint' is the new 'b'
        else
          a = midpoint;
      }
      
      const double tolerance = 100. * eps_double * std::abs((a+b))/2.;//eps gives rel. tol.
      if (show_output)
          std::cout << "\ttolerance = " << tolerance << std::endl;
      if (std::abs((b-a))<tolerance)
        converged = true; // Check if the half difference is below tolerance
      
    } while (!converged); // The looping stops if convergence has been reached.
    
    return (b+a)/2.; // Take mu to be the average of a and b.
}
