// kspace.cc
/* Description */

#include <iostream>
#include <cmath> // For constant M_PI
#include "alloc.h" // For allocation of multidimensional arrays
#include "array_init.h" // Initialization of arrays
#include "kspace.h" // Include header file for consistency check

// Internal subroutines


// Constructor implementation
kspace_t::kspace_t(const double* b1, const double* b2, const int b1_pts, const int b2_pts, 
                   const int bands_num, const bool with_output, const bool with_evecs)
    :bands_num_(bands_num), 
    b1_(b1), 
    b2_(b2), 
    b1_pts_(b1_pts), 
    b2_pts_(b2_pts),
    kx_grid(new double [b1_pts*b2_pts]), 
    ky_grid(new double [b1_pts*b2_pts]), 
    energies(Alloc2D_d(b1_pts*b2_pts, bands_num)),
    with_output_(with_output), 
    with_evecs_(with_evecs)
{
    if (with_output) std::cout << "kspace_t instance created.\n";
    // Assign MK momentum values
    MonkhorstPack_assign();
    
    // The energies grid is initialized to zero
    ValInitArray(b1_pts*b2_pts*bands_num, &(energies[0][0]), 0.);
    
    if (with_evecs) // It is convenient for evecs to be a 3D array.
    {
      evecs = Alloc3D_z(b1_pts*b2_pts, bands_num, bands_num);
      ValInitArray(b1_pts*b2_pts*bands_num*bands_num, &(evecs[0][0][0]), {0.,0.});
    }
}

// Destructor implementation
kspace_t::~kspace_t()
{
    delete [] kx_grid;
    delete [] ky_grid;
    Dealloc2D(energies);
    if (with_evecs_)
      Dealloc3D(evecs, b1_pts_*b2_pts_);
    if (with_output_) std::cout << "kspace_t instance deleted.\n";
}

int kspace_t::k_i(const int k1_ind, const int k2_ind) const
{
    // Return the global momentum index from the indices along each primitive lattice vector.
    return k2_ind + 
           b2_pts_*k1_ind;
}

void kspace_t::MonkhorstPack_assign()
{
    /* Assigns the MK momentum grid points to a 1D BZ given lattice cell length 'a'. We 
    use the notation of the MK article for simplicity. Choose num_pts even to avoid the 
    origin. 
    const int q = num_pts;
    for (int r=0; r<q; ++r)
        k_grid[r] = (M_PI/a)* (double)(2*r-q+1)/(double)(q);*/
    
    
    for (int i=0; i<b1_pts_; ++i)
      for (int j=0; j<b2_pts_; ++j)
      {
        kx_grid[k_i(i,j)] = (double)(2*i-b1_pts_+1)/(double)(2*b1_pts_)*b1_[0] + (double)(2*j-b2_pts_+1)/(double)(2*b2_pts_)*b2_[0];
        ky_grid[k_i(i,j)] = (double)(2*i-b1_pts_+1)/(double)(2*b1_pts_)*b1_[1] + (double)(2*j-b2_pts_+1)/(double)(2*b2_pts_)*b2_[1];
      }
}