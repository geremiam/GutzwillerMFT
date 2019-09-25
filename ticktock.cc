// ticktock.cc
//
// Supplies implementation of the class TickTock to measure wall-time
//
// Ramses van Zon, 2015-2016
//
// Modified by Geremia Massarelli in Jan 2018.

#include "ticktock.h"
#include <iostream>
#include <iomanip>

// checking for C++11 support, which has a standard timing library called chrono
#if __cplusplus > 199711L
  #include <chrono>
  // define a reference point:
  const static std::chrono::steady_clock::time_point reference_time = std::chrono::steady_clock::now();
#else
  // if not c++11, fall back on c calls from sys/time.h
  #include <sys/time.h>
#endif

static void print(const char* prefix, const double dt)
{
    // helper function to print elapsed time	
    if (prefix!=0)
        std::cout << prefix << "\t" << std::setprecision(4) << dt << " sec\n";
    else
        std::cout << "    \t" << std::setprecision(4) << dt << " sec\n";
}

static double elapsed_time()
{
  #if __cplusplus > 199711L
    std::chrono::steady_clock::time_point this_time = std::chrono::steady_clock::now();
    double dt = std::chrono::duration<double>(this_time - reference_time).count();
    return dt;
  #else
    struct timeval t;
    gettimeofday(&t, 0);
    return t.tv_sec + 0.000001 * t.tv_usec;
  #endif
}

void TickTock::tick()
{
    // start measuring time
    tick_ = elapsed_time();
}

void TickTock::tock(const char* prefix) const
{
    // done measuring time; print elapsed seconds
    print(prefix, silent_tock()); 
}

double TickTock::silent_tock() const
{
    // done measuring time; give elapsed seconds
    double tock = elapsed_time();
    double dt = tock-tick_;
    return dt;
}




Timer_t::Timer_t()
{
    // Constructor -- sets 'value_' to zero and 'running_' to false.
    value_ = 0.;
    running_ = false;
}

void Timer_t::start()
{
    if (!running_)
    {
        stopwatch.tick();
        running_ = true;
    }
    else
        std::cout << "Timer already running!" << std::endl;
}
void Timer_t::stop()
{
    if (running_)
    {
        value_ += stopwatch.silent_tock();
        running_ = false;
    }
    else
        std::cout << "Timer already stopped!" << std::endl;
}
double Timer_t::read()
{
    double reading=-99.;
    if (!running_)
        reading = value_;
    else
        std::cout << "Please stop timer before reading." << std::endl;
    
    return reading;
}
void Timer_t::reset()
{
    if (!running_)
        value_ = 0.;
    else
        std::cout << "Please stop timer before resetting." << std::endl;
}
