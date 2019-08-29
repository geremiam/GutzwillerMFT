// kspace.cc
/* Description */

#include <iostream>
#include <complex> // Needed to import alloc_dealloc.h
#include <cmath> // For constant M_PI
#include "alloc.h" // For allocation of multidimensional arrays
#include "init_routines.h" // Initialization of arrays
#include "kspace.h" // Include header file for consistency check

// Internal subroutines
void MonkhorstPack(const double a, const int num_pts, double*const k_grid)
{
    /* Assigns the MK momentum grid points to a 1D BZ given lattice cell length 'a'. We 
    use the notation of the MK article for simplicity. Choose num_pts even to avoid the 
    origin. */
    const int q = num_pts;
    for (int r=0; r<q; ++r)
        k_grid[r] = (M_PI/a)* (double)(2*r-q+1)/(double)(q);
}

// Constructor implementation
kspace_t::kspace_t(const double a, const double b, const double c, const int ka_pts, const int kb_pts, const int kc_pts, 
                   const int bands_num, const bool with_output, const bool with_evecs)
    :a_(a), ka_pts_(ka_pts), ka_grid(new double [ka_pts]), 
    b_(b), kb_pts_(kb_pts), kb_grid(new double [kb_pts]), 
    c_(c), kc_pts_(kc_pts), kc_grid(new double [kc_pts]), 
    bands_num_(bands_num), energies(new double [ka_pts*kb_pts*kc_pts*bands_num]),
    with_output_(with_output), with_evecs_(with_evecs)
{
    if (with_output) std::cout << "kspace_t instance created.\n";
    MonkhorstPack(a, ka_pts, ka_grid); // Assign MK momentum values along each axis
    MonkhorstPack(b, kb_pts, kb_grid);
    MonkhorstPack(c, kc_pts, kc_grid);
    // The energies grid is initialized to zero
    ValInitArray(ka_pts*kb_pts*kc_pts*bands_num, energies, 0.);
    
    if (with_evecs) // It is convenient for evecs to be a 3D array.
    {
      evecs = Alloc3D_z(ka_pts*kb_pts*kc_pts, bands_num, bands_num);
      ValInitArray(ka_pts*kb_pts*kc_pts*bands_num*bands_num, &(evecs[0][0][0]), {0.,0.});
    }
}

// Destructor implementation
kspace_t::~kspace_t()
{
    delete [] ka_grid;
    delete [] kb_grid;
    delete [] kc_grid;
    delete [] energies;
    if (with_evecs_)
      Dealloc3D(evecs, ka_pts_*kb_pts_*kc_pts_);
    if (with_output_) std::cout << "kspace_t instance deleted.\n";
}

int kspace_t::index(const int ka_ind, const int kb_ind, const int kc_ind, const int band_int) const
{
    /* Given indices for the momenta and bands, returns the corresponding index to be 
    used in the array "energies". */
    return band_int + 
           bands_num_*kc_ind + 
           bands_num_*kc_pts_*kb_ind + 
           bands_num_*kc_pts_*kb_pts_*ka_ind;
}

int kspace_t::k_ind(const int ka_ind, const int kb_ind, const int kc_ind) const
{
    // Return the global momentum index from the indices along each axis. To be used 
    // with "evecs".
    return kc_ind + 
           kc_pts_*kb_ind + 
           kc_pts_*kb_pts_*ka_ind;
}
