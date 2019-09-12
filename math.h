// math_routines.h
/*
Module containing various math functions or algorithms used to be used in other 
applications.
*/
#ifndef MATH_ROUTINES_H
#define MATH_ROUTINES_H

#include <complex> // Includes complex numbers
#include <limits> // To get machine precision

// Natural log of the largest double
const double log_max_double = std::log(std::numeric_limits<double>::max());
const double eps_double = std::numeric_limits<double>::epsilon(); // Used as rel. tol.


double log_1p_exp(const double x);

double nF0(double const energy);
float  nF0(float  const energy);
double nF(const double T, const double energy);
float  nF(const float  T, const float  energy);

double bisec(const double a_in, const double b_in, double foo(const double), const bool show_output=false);

#endif
