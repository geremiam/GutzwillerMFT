// kspace.h
/* Class for holding Monkhorst-Pack grid in momentum space.
Note that the object file must be linked to the alloc.o and array_init.o object files. */
#ifndef KSPACE_H
#define KSPACE_H

#include <complex>

class kspace_t
{
private:
    // Private copy constructor (prohibits copy creation)
    kspace_t(const kspace_t&);
    // Private assignment operator (prohibits assignment)
    const kspace_t& operator=(const kspace_t&);
    
    // Number of points along every reciprocal lattice vector. Choose even values to avoid the origin.
    const int b1_pts_;
    const int b2_pts_;
    
    // Reciprocal lattice vectors
    double b1_[2]={0.};
    double b2_[2]={0.};
    
    const int bands_num_; // Number of bands
    
    const bool show_output_; // Whether or not to silence diagnostic output
    const bool reserve_evals_; // Whether or not to reserve memory for evals
    const bool reserve_evecs_; // Whether or not to reserve memory for evecs
    
    void MonkhorstPack_assign(); // Assign MP momenta to the grids
    
public:
    /* The energies are held in a 1D array with the understanding that from slowest 
    varying to fastest varying, the indices are those for ka, kb, kc, and band.*/
    /* All momentum indices are collapsed into one; the function k_i is used to get it. */
    double*const kx_grid=NULL; // Only index is momentum
    double*const ky_grid=NULL; // Only index is momentum
    double*const* energies=NULL; // First index is momentum, second is band
    std::complex<double>*const*const* evecs=NULL; // Firs index is momentum; second and third are matrix indices
    
    // Constructor declaration
    kspace_t(const double* b1, const double* b2, const int b1_pts, const int b2_pts, 
             const int bands_num=0, const bool reserve_evals=false, const bool reserve_evecs=false, 
             const bool show_output=false);
    ~kspace_t(); // Destructor declaration
    
    // Return the global momentum index from the indices along each axis
    int k_i(const int k1_ind, const int k2_ind) const;
};

#endif
