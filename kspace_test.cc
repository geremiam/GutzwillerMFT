// kspace_test.cc
#include <iostream>
#include "kspace.h"

void kspace_t_test()
{
    const double a=1., b=1., c=1.;
    const int ka_pts=2, kb_pts=4, kc_pts = 6, bands_num=2;
    const bool with_output=false;
    const bool with_evecs=true;
    
    kspace_t Inst1(a, b, c, ka_pts, kb_pts, kc_pts, bands_num, with_output, with_evecs);
    
    /*
    std::cout << "Inst1.ka_pts_ = " << Inst1.ka_pts_ << std::endl;
    std::cout << "Inst1.kb_pts_ = " << Inst1.kb_pts_ << std::endl;
    std::cout << "Inst1.kc_pts_ = " << Inst1.kc_pts_ << std::endl;
    std::cout << "Inst1.bands_num_ = " << Inst1.bands_num_ << std::endl;
    */
    
    std::cout << "Inst1.ka_grid = " << std::endl;
    for (int i=0; i< ka_pts; ++i)
        std::cout << Inst1.ka_grid[i] << " ";
    std::cout << std::endl;
    
    std::cout << "Inst1.kb_grid = " << std::endl;
    for (int i=0; i< kb_pts; ++i)
        std::cout << Inst1.kb_grid[i] << " ";
    std::cout << std::endl;
    
    std::cout << "Inst1.kc_grid = " << std::endl;
    for (int i=0; i< kc_pts; ++i)
        std::cout << Inst1.kc_grid[i] << " ";
    std::cout << std::endl;
    
    // Assign values to energies
    for (int i=0 ; i<ka_pts*kb_pts*kc_pts*bands_num; ++i)
      Inst1.energies[i] = (double)(i);
    // Assign values to evecs
    for (int i=0; i<ka_pts*kb_pts*kc_pts*bands_num*bands_num; ++i)
      Inst1.evecs[0][0][i] = {(double)(i), 0.1*(double)(i)};
    
    // Print them out using the 'index()' function
    std::cout << "Inst1.energies = " << std::endl;
    for (int i=0; i<ka_pts; ++i)
    {
      for (int j=0; j<kb_pts; ++j)
      {
        for (int k=0; k<kc_pts; ++k)
        {
          for (int l=0; l<bands_num; ++l)
            std::cout << Inst1.energies[Inst1.index(i,j,k,l)] << " ";
            //std::cout << Inst1.evecs[Inst1.k_ind(i,j,k)][l][0] << " ";
          std::cout << "\t";
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
