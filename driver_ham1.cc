// driver_ham1.cc
/* Mean-field theory of antiferromagnetism in the Haldane bilayer, where layer II also 
has a Hubbard interaction. This driver defines model-specific parameters and functions, 
and performs the iteration until self-consistency is achieved. The data is saved as a 
NetCDF dataset. */
#include <iostream>
#include <sstream> // For stringstreams
#include <complex> // For complex numbers
#include <cmath> // For many math functions
#include <string>
#include <cassert>
#include "alloc.h"
#include "array_init.h" // Array initialization
#include "ncio.h" // Class for creating simple NetCDF datasets
#include "IO.h" // For PrintMatrix()
#include "ham1.h" // Source code for ham1
using std::complex;
using std::string;

// Class that defines the parameter space for this Hamiltonian

class pspacebase_t {
  private:
    // Private copy constructor (prohibits copy creation)
    pspacebase_t(const pspacebase_t&);
    // Private assignment operator (prohibits assignment)
    const pspacebase_t& operator=(const pspacebase_t&);
    
  public:
    // Parameter defining the parameter space
    
    // The order of array components must match.
    const size_t dims_num_; // Number of coordinates; to be given in constructor.
    string* dim_names_; // Name of every coordinate; to be given in constructor.
    size_t* dim_lengths_; // Number of points for every coordinate, to be given in constructor.
    double* dim_ranges_; // Range of the coordinates, given in the form [min0,max0,min1,max1,...] in the constructor.
    
    int parspace_pts; // To contain the number of points in parameter space (product of all the dim_lengths_).
    
    // Arrays of coordinate variables
    double** coord_arrays; // Each element is a pointer to an array of doubles.
    
    // Variables to be stored on the parameter space. We store them in 1D arrays, and use 
    // the method idx() to index them.
    // Mean fields
    double* chi_s_grid;
    double* chi_d_grid;
    complex<double>* Delta_s_grid;
    complex<double>* Delta_d_grid;
    // Other variables
    complex<double>* DeltaSC_s_grid;
    complex<double>* DeltaSC_d_grid;
    double* optweight_xx_grid;
    double* optweight_yy_grid;
    double* tp_grid;
    
    int   * loops_grid; // holds the number of loops done at each point
    double* energy_grid; // Holds the MF energy for later comparison
    double* mu_grid; // Holds the chemical potential for the converged system
    
    // Constructor declaration
    pspacebase_t(const size_t dims_num, const char*const dim_names[], const size_t*const dim_lengths, const double*const dim_ranges)
    :dims_num_(dims_num)
    {
        // Arrays that hold information on dimensions
        dim_names_   = new string   [dims_num_];
        dim_lengths_ = new size_t   [dims_num_];
        dim_ranges_  = new double [2*dims_num_];
        coord_arrays = new double*  [dims_num_];
        
        parspace_pts = 1; // Initialize to one to begin with
        
        // Loop over all dimensions
        for (int i=0; i<dims_num_; ++i)
        {
          if (dim_lengths[i]<1) std::cout << "*** WARNING *** ALL DIMENSION LENGTHS SHOULD BE GREATER THAN ZERO\n"; 
          // Copy values from input arrays into class attributes
          dim_names_[i]      = dim_names[i];
          dim_lengths_[i]    = dim_lengths[i];
          dim_ranges_[2*i]   = dim_ranges[2*i];
          dim_ranges_[2*i+1] = dim_ranges[2*i+1];
          
          parspace_pts *= dim_lengths_[i]; // Find total size of parspace
          
          // Allocate and initialize arrays of coordinates based on given values
          coord_arrays[i] = new double [dim_lengths_[i]]; // Allocate array
          if (dim_lengths_[i]>1) // Initialize multi-component arrays
          {
            const bool endpoint = true; // (not a particularly important choice)
            LinInitArray(dim_ranges_[2*i], dim_ranges_[2*i+1], dim_lengths_[i], coord_arrays[i], endpoint);
          }
          else // Initialize single-component arrays
            coord_arrays[i][0] = dim_ranges_[2*i];
        }
        
        // Allocation of memory for the variables, which span the whole coordinate space
        chi_s_grid        = new double          [parspace_pts];
        chi_d_grid        = new double          [parspace_pts];
        Delta_s_grid      = new complex<double> [parspace_pts];
        Delta_d_grid      = new complex<double> [parspace_pts];
        
        DeltaSC_s_grid    = new complex<double> [parspace_pts];
        DeltaSC_d_grid    = new complex<double> [parspace_pts];
        optweight_xx_grid = new double          [parspace_pts];
        optweight_yy_grid = new double          [parspace_pts];
        tp_grid           = new double          [parspace_pts];
        loops_grid        = new int             [parspace_pts];
        energy_grid       = new double          [parspace_pts];
        mu_grid           = new double          [parspace_pts];
        
        // Initialize variables to unlikely values.
        ValInitArray(parspace_pts, chi_s_grid,   -99.);
        ValInitArray(parspace_pts, chi_d_grid,   -99.);
        ValInitArray(parspace_pts, Delta_s_grid, -99.);
        ValInitArray(parspace_pts, Delta_d_grid, -99.);
        
        ValInitArray(parspace_pts, DeltaSC_s_grid, -99.);
        ValInitArray(parspace_pts, DeltaSC_d_grid, -99.);
        ValInitArray(parspace_pts, optweight_xx_grid, -99.);
        ValInitArray(parspace_pts, optweight_yy_grid, -99.);
        ValInitArray(parspace_pts, tp_grid, -99.);
        ValInitArray(parspace_pts, loops_grid,     -1);
        ValInitArray(parspace_pts, energy_grid,  -99.);
        ValInitArray(parspace_pts, mu_grid,       -9.);
        
        std::cout << "pspacebase_t instance created.\n";
    }
    // Destructor declaration
    ~pspacebase_t()
    {
        delete [] chi_s_grid; // MF variables
        delete [] chi_d_grid;
        delete [] Delta_s_grid;
        delete [] Delta_d_grid;
        
        delete [] DeltaSC_s_grid;
        delete [] DeltaSC_d_grid;
        delete [] optweight_xx_grid;
        delete [] optweight_yy_grid;
        delete [] tp_grid;
        delete [] loops_grid;
        delete [] energy_grid;
        delete [] mu_grid;
        
        for (int i=0; i<dims_num_; ++i)
          delete [] coord_arrays[i];
        
        delete [] dim_names_;
        delete [] dim_lengths_;
        delete [] dim_ranges_;
        delete [] coord_arrays; // Make sure the inner pointers are deleted before this one
        
        std::cout << "pspacebase_t instance deleted.\n";
    }
    
    virtual void parvals(const int i, double*const t_pointer, double*const tp_pointer, double*const J_pointer, double*const x_pointer, const bool with_output=false)
    {
        std::cout << "WARNING: Called pspacebase_t::parvals()\n";
    }
    
    void SaveData(const string GlobalAttr, const string path)
    {
        // Method for saving the data of this class. This method uses the class from the 
        // module nc_IO that creates a simple NetCDF class and allows writing of variables.
        // We define parameters required to create the dataset. Don't forget to adjust these depending on the parameter space defined above. 
        // Important: Note that the order of the variables must be kept consistent.
        
        const size_t vars_num = 11; // Variables to be saved other than coord variables
        string var_names [vars_num] = {"chi_s","chi_d","Delta_s","Delta_d","DeltaSC_s","DeltaSC_d","optweight_xx","optweight_yy", "tp", "mu","energy"}; // List for the variable names
        bool var_complex [vars_num] = {  false,  false,     true,     true,       true,       true,         false,         false,false,false,false}; // List for indicating whether vars are complex
        
        // Constructor for the dataset class creates a dataset
        const bool prune = true; // Activates pruning of length-1 dimensions.
        newDS_t newDS(dims_num_, dim_names_, dim_lengths_, vars_num, var_names, var_complex, GlobalAttr, path, prune);
        
        // Define coordinate variables. Order matters.
        for (int i=0; i<dims_num_; ++i)
          newDS.DefCoordVar(i, dim_names_[i]);
        
        newDS.EndDef(); // Exit definition mode
        
        // Write coordinate variables
        for (int i=0; i<dims_num_; ++i)
          newDS.WriteCoordVar(i, coord_arrays[i]);
        
        // List for holding the pointers to the vars
        double* vars [vars_num] = {chi_s_grid, chi_d_grid, 
                                   reinterpret_cast<double*const>(Delta_s_grid), 
                                   reinterpret_cast<double*const>(Delta_d_grid),
                                   reinterpret_cast<double*const>(DeltaSC_s_grid), 
                                   reinterpret_cast<double*const>(DeltaSC_d_grid),
                                   optweight_xx_grid, optweight_yy_grid, tp_grid,
                                   mu_grid,
                                   energy_grid};
        newDS.WriteVars(vars); // Write the variables
        
        newDS.WriteLoops(loops_grid); // Write loops variable
    }
    
    int pstudy(const MFs_t MFs_initial, FPparams_t FPparams, const int k1_pts = 298, const int k2_pts = 298, const bool with_output = false)
    {
        // Choose twice a prime number for momentum space grid resolution
        // This routine performs the mean-field iterative search at every point in the 
        // parameter space defined above. Show output for diagnostics if with_output.
        int numfails = 0; // Tracks number of points which failed to converge after loops_lim
        
        string GlobalAttr; // String in which to save the global attributes.
        
        // Must declare 'const' variables as firstprivate for compatibility with new gcc 
        // versions. See https://gcc.gnu.org/gcc-9/porting_to.html#ompdatasharing.
        #pragma omp parallel default(none) shared(FPparams,GlobalAttr,std::cout) firstprivate(with_output,k1_pts,k2_pts,MFs_initial) reduction(+:numfails)
        {
        // Declare and construct an instance of ham1_t
        ham1_t ham1(FPparams, MFs_initial, k1_pts, k2_pts);
    
        // Set parameters to start values. Useful to put in attributes. Leave temperature at default value.
        parvals(0, &ham1.t_, &ham1.tp_, &ham1.J_, &ham1.x_, with_output);
        
        #pragma omp single
        {
          GlobalAttr = ham1.GetAttributes(); // assign attributes to GlobalAttr
          
          if (with_output)
            std::cout << "\n\nFPparams:\n" << FPparams << "\n"
                      << "*********************************************************" << "\n";
        }
    
        // Loop over values of the parameter space
        #pragma omp for schedule(dynamic,1)
        for (int i=0; i<parspace_pts; ++i)
           {
               // Adjust phase space parameters
               parvals(i, &ham1.t_, &ham1.tp_, &ham1.J_, &ham1.x_, with_output);
               
               int loops=0; // Will be assigned the number of loops performed
               double mu=0.;// Will be assigned the chemical potential of the converged Ham.
               double energy=0.;
               complex<double> DeltaSC_s=0.;
               complex<double> DeltaSC_d=0.;
               double optweight_xx=0.;
               double optweight_yy=0.;
               const bool fail = ham1.FixedPoint(with_output, &loops, &mu, &energy, &DeltaSC_s, &DeltaSC_d, &optweight_xx, &optweight_yy);
               
               if (fail) // Print current params
               {
                   std::cout << "\tWARNING: failure to converge after limit reached at point " << i << "of " << parspace_pts
                             <<  "\tt = " << ham1.t_
                             << "\ttp = " << ham1.tp_  
                             <<  "\tJ = " << ham1.J_ 
                             <<  "\tx = " << ham1.x_ << "\n";
                   ++numfails;
               }
               
               chi_s_grid  [i] = ham1.MFs_.chi_s;
               chi_d_grid  [i] = ham1.MFs_.chi_d;
               Delta_s_grid[i] = ham1.MFs_.Delta_s;
               Delta_d_grid[i] = ham1.MFs_.Delta_d;
               
               DeltaSC_s_grid   [i] = DeltaSC_s;
               DeltaSC_d_grid   [i] = DeltaSC_d;
               optweight_xx_grid[i] = optweight_xx;
               optweight_yy_grid[i] = optweight_yy;
               tp_grid          [i] = ham1.tp_;
               energy_grid      [i] = energy;
               mu_grid          [i] = mu;
               loops_grid       [i] = loops; // Save the number of loops to array.
               
               if (with_output) std::cout << std::endl;
           }
        }
    
        // We save to a NetCDF dataset using the class defined in nc_IO. Call saving method.
        SaveData(GlobalAttr, "data/ham1/"); // Include final '/' in path for saving
    
        // Print out numfails
        std::cout << "\n\tNUMBER OF FAILURES TO CONVERGE: " << numfails << "\n\n";
    
        return numfails;
    }
};

class pspace_t:public pspacebase_t {
  private:
    // Private copy constructor (prohibits copy creation)
    pspace_t(const pspace_t&);
    // Private assignment operator (prohibits assignment)
    const pspace_t& operator=(const pspace_t&);
    
  public:
    // Parameter defining the parameter space
    
    // The order of array components must match.
    // Order of coordinates: t, tp, J, x
    static const size_t dims_num = 4; // Number of coordinates
    static constexpr const char*const dim_names[dims_num] = {"t_", "tp_", "J_", "x_"};
    
    // Constructor declaration
    pspace_t(const size_t*const dim_lengths, const double*const dim_ranges)
    :pspacebase_t(dims_num, dim_names, dim_lengths, dim_ranges)
    {
        std::cout << "pspace_t instance created.\n";
    }
    // Destructor declaration
    ~pspace_t()
    {
        std::cout << "pspace_t instance deleted.\n";
    }
    
    void parvals(const int i, double*const t_pointer, double*const tp_pointer, double*const J_pointer, double*const x_pointer, const bool with_output=false)
    {
        // Returns values of the parameters based on the index i. Indexing is in the order
        // of the variables given in dim_names_, etc.
        
        // Recall the order: t_, tp_, J_, x_
        
        // Define a reference for convenience
        const size_t*const& len = dim_lengths_;
        
        if ( (i >= parspace_pts) || (i < 0) )
          std::cout << "WARNING from method parvals(): i is out of range.\n";
        
        const int i3 = i % len[3];
        const int i2 = ( (i - i3) / len[3] ) % len[2];
        const int i1 = ( (i - i3 - i2*len[3]) / (len[2]*len[3]) ) % len[1];
        const int i0 = ( (i - i3 - i2*len[3] - i1*len[2]*len[3]) / (len[1]*len[2]*len[3]) ) % len[0];
        
        assert( i == i3 + 
                     i2 * len[3] + 
                     i1 * len[3] * len[2] + 
                     i0 * len[3] * len[2] * len[1]);
        
        // Assign values
        *t_pointer  = coord_arrays[0][i0];
        *tp_pointer = coord_arrays[1][i1];
        *J_pointer  = coord_arrays[2][i2];
        *x_pointer  = coord_arrays[3][i3];
        
        if (with_output) // Print current params
        {
            std::cout << "\n\n";
            if (len[0]!=1) std::cout <<  "t = " << *t_pointer  << "\t";
            if (len[1]!=1) std::cout << "tp = " << *tp_pointer << "\t";
            if (len[2]!=1) std::cout <<  "J = " << *J_pointer  << "\t";
            if (len[3]!=1) std::cout <<  "x = " << *x_pointer  << "\t";
            std::cout << "\n";
        }
    }
    
};

class pspace_time_t:public pspacebase_t {
  private:
    // Private copy constructor (prohibits copy creation)
    pspace_time_t(const pspace_time_t&);
    // Private assignment operator (prohibits assignment)
    const pspace_time_t& operator=(const pspace_time_t&);
    
  public:
    // Parameter defining the parameter space
    
    // The order of array components must match.
    // Order of coordinates: t, tp_avg, tp_amp, time, J, x 
    // coord_arrays for time is not used (only index is used)
    static const size_t dims_num = 6; // Number of coordinates
    static constexpr const char*const dim_names[dims_num] = {"t_", "tp_avg", "tp_amp", "time", "J_", "x_"};
    
    // Constructor declaration
    pspace_time_t(const size_t*const dim_lengths, const double*const dim_ranges)
    :pspacebase_t(dims_num, dim_names, dim_lengths, dim_ranges)
    {
        std::cout << "pspace_time_t instance created.\n";
    }
    // Destructor declaration
    ~pspace_time_t()
    {
        std::cout << "pspace_time_t instance deleted.\n";
    }
    
    void parvals(const int i, double*const t_pointer, double*const tp_pointer, double*const J_pointer, double*const x_pointer, const bool with_output=false)
    {
        // Returns values of the parameters based on the index i. Indexing is in the order
        // of the variables given in dim_names_, etc.
        
        // Recall the order: t_, tp_avg, tp_amp, time, J_, x_
        
        // Define a reference for convenience
        const size_t*const& len = dim_lengths_;
        
        if ( (i >= parspace_pts) || (i < 0) )
          std::cout << "WARNING from method parvals(): i is out of range.\n";
        
        const int i5 = i % len[5];
        const int i4 = ( (i - i5) / len[5] ) % len[4];
        const int i3 = ( (i - i5 - i4*len[5]) / (len[4]*len[5]) ) % len[3];
        const int i2 = ( (i - i5 - i4*len[5] - i3*len[4]*len[5]) / (len[3]*len[4]*len[5]) ) % len[2];
        const int i1 = ( (i - i5 - i4*len[5] - i3*len[4]*len[5] - i2*len[3]*len[4]*len[5]) / (len[2]*len[3]*len[4]*len[5]) ) % len[1];
        const int i0 = ( (i - i5 - i4*len[5] - i3*len[4]*len[5] - i2*len[3]*len[4]*len[5] - i1*len[2]*len[3]*len[4]*len[5]) / (len[1]*len[2]*len[3]*len[4]*len[5]) ) % len[0];
        
        assert( i == i5 + 
                     i4 * len[5] + 
                     i3 * len[4] * len[5] + 
                     i2 * len[3] * len[4] * len[5] + 
                     i1 * len[2] * len[3] * len[4] * len[5] + 
                     i0 * len[1] * len[2] * len[3] * len[4] * len[5]);
        
        // Assign values
        *t_pointer  = coord_arrays[0][i0];
        //             (__ tp_avg __)       (__ tp_amp ___)                    (_ time _) (period)
        *tp_pointer = coord_arrays[1][i1] + coord_arrays[2][i2] * sin(2.*M_PI*(double)(i3)/len[3]);
        // Note that there is no double counting because i3 goes from 0 to len[3]-1.
        // coord_arrays for time is not used; instead, index is used.
        *J_pointer  = coord_arrays[4][i4];
        *x_pointer  = coord_arrays[5][i5];
        
        if (with_output) // Print current params
        {
            std::cout <<  "\n\nt_idx = " << i0 << " of " << len[0]-1 << "    ";
            std::cout << "tp_avg_idx = " << i1 << " of " << len[1]-1 << "    ";
            std::cout << "tp_amp_idx = " << i2 << " of " << len[2]-1 << "    ";
            std::cout <<   "time_idx = " << i3 << " of " << len[3]-1 << "    ";
            std::cout <<      "J_idx = " << i4 << " of " << len[4]-1 << "    ";
            std::cout <<      "x_idx = " << i5 << " of " << len[5]-1 << "\n";
            
            std::cout <<  "t = " << *t_pointer  << "    ";
            std::cout << "tp = " << *tp_pointer << "    ";
            std::cout <<  "J = " << *J_pointer  << "    ";
            std::cout <<  "x = " << *x_pointer  << "\n\n";
        }
    }
    
};

// ***************************************************************************************
int pspace_studyA(const bool show_output)
{
    // Order of coordinates: t, tp, J, x
    const size_t dim_lengths[4]  = {1, 1, 1, 61};
    const double dim_ranges[2*4] = {1., 1.,    0., 0.,    1./5., 1./5.,    0.01, 0.31};
    
    pspace_t pspace(dim_lengths, dim_ranges);
    
    MFs_t MFs_initial(0.1, 0., {0.,0.}, {0.1,0.});
    
    const double tol = 1.e-6;
    const int loops_lim = 800;
    const int mixing_vals_len = 7;
    const int   counter_vals [mixing_vals_len] = {100, 200, 300, 400, 500, 600, 700};
    const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3};
    FPparams_t FPparams(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    
    // Choose twice a prime number for momentum space grid resolution
    const int k1_pts = 298;
    const int k2_pts = 298;
    
    const int info = pspace.pstudy(MFs_initial, FPparams, k1_pts, k2_pts, show_output);
    return info;
}

int pspace_studyB(const bool show_output)
{
    // Order of coordinates: t, tp, J, x
    const size_t dim_lengths[4]  = {1, 101, 1, 31}; // {1, 26, 1, 16};
    const double dim_ranges[2*4] = {1., 1.,    -0.5, 0.5,    1./5., 1./5.,    0.01, 0.31};
    
    pspace_t pspace(dim_lengths, dim_ranges);
    
    MFs_t MFs_initial(0.1, 0., {0.,0.}, {0.1,0.});
    
    const double tol = 1.e-6;
    const int loops_lim = 800;
    const int mixing_vals_len = 7;
    const int   counter_vals [mixing_vals_len] = {100, 200, 300, 400, 500, 600, 700};
    const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3};
    FPparams_t FPparams(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    
    // Choose twice a prime number for momentum space grid resolution
    const int k1_pts = 298;
    const int k2_pts = 298;
    
    const int info = pspace.pstudy(MFs_initial, FPparams, k1_pts, k2_pts, show_output);
    return info;
}
// ***************************************************************************************
int pspace_time_studyA(const bool show_output)
{
    // Order of coordinates:        t,      tp_avg,   tp_amp,   time,   J,            x 
    const size_t dim_lengths[6]  = {1,      1,        5,        100,    1,            16};
    const double dim_ranges[2*6] = {1., 1., -.3, -.3, .01, .05, 0., 1., 1./5., 1./5., 0.01, 0.31};
    
    pspace_time_t pspace(dim_lengths, dim_ranges);
    
    MFs_t MFs_initial(0.1, 0., {0.,0.}, {0.1,0.});
    
    const double tol = 1.e-6;
    const int loops_lim = 800;
    const int mixing_vals_len = 7;
    const int   counter_vals [mixing_vals_len] = {100, 200, 300, 400, 500, 600, 700};
    const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3};
    FPparams_t FPparams(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    
    // Choose twice a prime number for momentum space grid resolution
    const int k1_pts = 298;
    const int k2_pts = 298;
    
    const int info = pspace.pstudy(MFs_initial, FPparams, k1_pts, k2_pts, show_output);
    return info;
}
// ***************************************************************************************


int main(int argc, char* argv[])
{
    const int info = pspace_time_studyA(false);
    return info;
}
