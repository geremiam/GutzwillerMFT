// ham1_test.cc
#include <iostream>
#include <iomanip> // For the function std::setprecision()
#include <complex> // For complex numbers
#include <cmath> // For the constant M_PI
#include "IO.h" // For printing matrices
#include "alloc.h"
#include "ham1.h"

void Assign_ham_test()
{
    std::complex<double>** matrix = Alloc2D_z(8, 8);
    std::cout << "\nmatrix = " << std::endl;
    PrintMatrix(8,8,matrix, std::cout);
    std::cout << std::endl;
    
    const int ka_pts = 10;
    const int kb_pts = 10;
    const int kc_pts = 10;
    
    ham1_t ham1(ka_pts, kb_pts, kc_pts);
    
    std::cout << "num_states = " << ham1.num_states << std::endl;
    std::cout << "rho_ = " << ham1.rho_ << std::endl;
    std::cout << "filled_states = " << ham1.filled_states << std::endl;
    
    const double ka=0.1*M_PI, kb=0.1*M_PI, kc=0.1*M_PI;
    ham1.Assign_ham(ka, kb, kc, matrix);
    std::cout << "\nAfter assignment, matrix = " << std::endl;
    PrintMatrix(8,8,matrix, std::cout);
    std::cout << std::endl;
    
    /*
    ham1.Assign_ham(0.,0.);
    PrintMatrix(8,8,ham1.ham_array,std::cout);
    
    std::cout << std::endl;
    
    ham1.U = 0.5;
    ham1.Assign_ham(0.,0.);
    PrintMatrix(8,8,ham1.ham_array,std::cout);
    */
    Dealloc2D(matrix);
}

void test_ComputeMFs()
{
    const int ka_pts = 62;
    const int kb_pts = 62;
    const int kc_pts = 62;
    
    ham1_t ham1(ka_pts, kb_pts, kc_pts);
    ham1.set_nonzerotemp(1.e-1);
    //ham1.set_zerotemp();
    
    /*ham1.assign_rho(0.7);
    ham1.V1_ = 40.;
    ham1.V1p_ = 30.;
    ham1.V2_ = 6.;
    ham1.V3_ = 8.;
    
    ham1.rho_s_[0] = 0.2;    ham1.rho_a_[0] = 0.1;
    ham1.rho_s_[1] = 0.3;    ham1.rho_a_[1] = 0.1;
    ham1.rho_s_[2] = 0.0;    ham1.rho_a_[2] = 0.0;
    ham1.rho_s_[3] = 0.0;    ham1.rho_a_[3] = 0.0;*/
    
    std::cout << "num_states: " << ham1.num_states << std::endl;
    std::cout << "filled_states: " << ham1.filled_states << std::endl;
    
    double rho_s_out [4] = {0.};
    double rho_a_out [4] = {0.};
    
    ham1.ComputeMFs_old(rho_s_out, rho_a_out);
    for (int i=0; i<4; ++i) {
      std::cout << "rho_s_out[" << i << "] = " << rho_s_out[i] << "\t"
                << "rho_a_out[" << i << "] = " << rho_a_out[i] << std::endl;
    }
    
    ham1.ComputeMFs(rho_s_out, rho_a_out);
    for (int i=0; i<4; ++i) {
      std::cout << "rho_s_out[" << i << "] = " << rho_s_out[i] << "\t"
                << "rho_a_out[" << i << "] = " << rho_a_out[i] << std::endl;
    }
}

void test_FixedPoint()
{
    const int ka_pts = 62;
    const int kb_pts = 62;
    const int kc_pts = 62;
    
    ham1_t ham1(ka_pts, kb_pts, kc_pts);
    ham1.set_zerotemp(); // In this case, calculates EF from sorting
    //ham1.set_nonzerotemp(1.e-2); // In this case, calculates mu from bisection method
    
    int num_loops;
    const bool with_output = true;
    
    ham1.assign_rho(0.7);
    ham1.V1_ = 40.;
    ham1.V1p_ = 30.;
    ham1.V2_ = 6.;
    ham1.V3_ = 8.;
    
    ham1.FixedPoint(&num_loops, with_output);
    
    std::cout << "num_loops = " << num_loops << std::endl;
}

int main()
{
    //Assign_ham_test();
    
    //test_ComputeMFs();
    
    test_FixedPoint();
    
    return 0;
}
