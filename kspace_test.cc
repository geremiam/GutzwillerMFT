// kspace_test.cc
#include <iostream>
#include "kspace.h"

void kspace_t_test()
{
    const int bands_num=2;
    const bool with_output=false;
    const bool with_evals=true;
    const bool with_evecs=true;
    
    const double b1[2] = {1., 0.};
    const double b2[2] = {0., 1.};
    
    // Addresses of arrays
    std::cout << "b1 = " << b1 << std::endl;
    std::cout << "b2 = " << b2 << std::endl;
    
    const int b1_pts = 3;
    const int b2_pts = 4;
    
    kspace_t Inst1(b1, b2, b1_pts, b2_pts, bands_num, with_evals, with_evecs, with_output);
    
    // Show parameters
//     std::cout << "Inst1.b1_ = " << Inst1.b1_ << std::endl;
//     std::cout << "Inst1.b2_ = " << Inst1.b2_ << std::endl;
//     std::cout << "Inst1.b1_pts_ = " << Inst1.b1_pts_ << std::endl;
//     std::cout << "Inst1.b2_pts_ = " << Inst1.b2_pts_ << std::endl;
//     std::cout << "Inst1.bands_num_ = " << Inst1.bands_num_ << std::endl;
    
    // Show momentum grids
    std::cout << "Inst1.kx_grid = " << std::endl;
    for (int i=0; i<b1_pts; ++i)
    {
      for (int j=0; j<b2_pts; ++j)
        std::cout << Inst1.kx_grid[Inst1.k_i(i,j)] << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
    
    std::cout << "Inst1.ky_grid = " << std::endl;
    for (int i=0; i<b1_pts; ++i)
    {
      for (int j=0; j<b2_pts; ++j)
        std::cout << Inst1.ky_grid[Inst1.k_i(i,j)] << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
    
    
    // Assign values to energies
    for (int i=0 ; i<b1_pts*b2_pts*bands_num; ++i)
      Inst1.energies[0][i] = (double)(i);
    // Assign values to evecs
    for (int i=0; i<b1_pts*b2_pts*bands_num*bands_num; ++i)
      Inst1.evecs[0][0][i] = {(double)(i), 0.1*(double)(i)};
    
    
    // Print them out using the 'index()' function
    std::cout << "Inst1.energies = " << std::endl;
    for (int i=0; i<b1_pts; ++i)
    {
      for (int j=0; j<b2_pts; ++j)
      {
          for (int l=0; l<bands_num; ++l)
          {
            std::cout << Inst1.energies[Inst1.k_i(i,j)][l] << " ";
            std::cout << Inst1.evecs[Inst1.k_i(i,j)][l][0] << "\t";
          }
          std::cout << "\n";
      }
      std::cout << "\n\n";
    }
}

int main()
{
    kspace_t_test();
    return 0;
}
