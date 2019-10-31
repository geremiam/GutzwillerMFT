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
#include "alloc.h"
#include "array_init.h" // Array initialization
#include "ncio.h" // Class for creating simple NetCDF datasets
#include "IO.h" // For PrintMatrix()
#include "ham1.h" // Source code for ham1
using std::complex;
using std::string;

// Class that defines the parameter space for this Hamiltonian

class pspace_t {
  private:
    // Private copy constructor (prohibits copy creation)
    pspace_t(const pspace_t&);
    // Private assignment operator (prohibits assignment)
    const pspace_t& operator=(const pspace_t&);
    
  public:
    // Parameter defining the parameter space
    
    // Fix the order: t_, tp_, J_, x_
    static const size_t dims_num = 4; // Number of coordinates
    const string dim_names[dims_num] = {"t_", "tp_", "J_", "x_"};
    
    size_t dim_lengths_  [dims_num]; // Number of points for every coordinate, to be given by user.
    double dim_ranges_ [2*dims_num]; // Range of the coordinates, given in the form [min0,max0,min1,max1,...]. To be given by user.
    
    int parspace_pts; // To contain the number of points in parameter space.
    
    // Arrays of coordinate variables
    double* coord_arrays [dims_num]; // Each element is a pointer to an array of doubles.
    
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
    
    int   * loops_grid; // holds the number of loops done at each point
    double* energy_grid; // Holds the MF energy for later comparison
    double* mu_grid; // Holds the chemical potential for the converged system
    
    // Constructor declaration
    pspace_t(const size_t*const dim_lengths, const double*const dim_ranges)
    {
        parspace_pts = 1; // Initialize to one to begin with
        
        for (int i=0; i<dims_num; ++i)
        {
          if (dim_lengths[i]<1) std::cout << "*** WARNING *** ALL DIMENSION LENGTHS SHOULD BE GREATER THAN ZERO\n"; 
          // Copy values from input arrays into class attributes
          dim_lengths_[i]    = dim_lengths[i];
          dim_ranges_[2*i]   = dim_ranges[2*i];
          dim_ranges_[2*i+1] = dim_ranges[2*i+1];
          
          parspace_pts *= dim_lengths_[i]; // Find total size of parspace
          
          // Initialize arrays of coordinates based on given values
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
        ValInitArray(parspace_pts, loops_grid,     -1);
        ValInitArray(parspace_pts, energy_grid,  -99.);
        ValInitArray(parspace_pts, mu_grid,       -9.);
        
        std::cout << "pspace_t instance created.\n";
    }
    // Destructor declaration
    ~pspace_t()
    {
        for (int i=0; i<dims_num; ++i)
          delete [] coord_arrays[i];
        
        delete [] chi_s_grid; // MF variables
        delete [] chi_d_grid;
        delete [] Delta_s_grid;
        delete [] Delta_d_grid;
        
        delete [] DeltaSC_s_grid;
        delete [] DeltaSC_d_grid;
        delete [] optweight_xx_grid;
        delete [] optweight_yy_grid;
        delete [] loops_grid;
        delete [] energy_grid;
        delete [] mu_grid;
        
        std::cout << "pspace_t instance deleted.\n";
    }
    
    int idx(const int t_idx, const int tp_idx, const int J_idx, const int x_idx)
    {
        // Returns the index of the var arrays corresponding to indices i, j of the coord
        // var arrays. Trivial in this case.
        return  x_idx + 
                J_idx * dim_lengths_[3] + 
               tp_idx * dim_lengths_[3] * dim_lengths_[2] + 
                t_idx * dim_lengths_[3] * dim_lengths_[2] * dim_lengths_[1];
    }
    
    void SaveData(const string GlobalAttr, const string path)
    {
        // Method for saving the data of this class. This method uses the class from the 
        // module nc_IO that creates a simple NetCDF class and allows writing of variables.
        // We define parameters required to create the dataset. Don't forget to adjust these depending on the parameter space defined above. 
        // Important: Note that the order of the variables must be kept consistent.
        
        const size_t vars_num = 10; // Variables to be saved other than coord variables
        string var_names [vars_num] = {"chi_s", "chi_d", "Delta_s", "Delta_d", "DeltaSC_s", "DeltaSC_d", "optweight_xx", "optweight_yy", "mu", "energy"}; // List for the variable names
        bool var_complex [vars_num] = { false,  false,   true,      true,      true,        true,        false,          false,          false, false}; // List for indicating whether vars are complex
        
        // Constructor for the dataset class creates a dataset
        const bool prune = true; // Activates pruning of length-1 dimensions.
        newDS_t newDS(dims_num, dim_names, dim_lengths_, vars_num, var_names, var_complex, GlobalAttr, path, prune);
        
        // Define coordinate variables. Order matters.
        for (int i=0; i<dims_num; ++i)
          newDS.DefCoordVar(i, dim_names[i]);
        
        newDS.EndDef(); // Exit definition mode
        
        // Write coordinate variables
        for (int i=0; i<dims_num; ++i)
          newDS.WriteCoordVar(i, coord_arrays[i]);
        
        // List for holding the pointers to the vars
        double* vars [vars_num] = {chi_s_grid, chi_d_grid, 
                                   reinterpret_cast<double*const>(Delta_s_grid), 
                                   reinterpret_cast<double*const>(Delta_d_grid),
                                   reinterpret_cast<double*const>(DeltaSC_s_grid), 
                                   reinterpret_cast<double*const>(DeltaSC_d_grid),
                                   optweight_xx_grid, optweight_yy_grid,
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
        #pragma omp parallel default(none) shared(FPparams,MFs_initial,GlobalAttr,std::cout) firstprivate(with_output) reduction(+:numfails)
        {
        // Declare and construct an instance of ham1_t
        ham1_t ham1(FPparams, MFs_initial, k1_pts, k2_pts);
    
        // Set parameters to start values. Needed in the attributes. Leave temperature at default value.
        ham1.t_  = dim_ranges_[2*0];
        ham1.tp_ = dim_ranges_[2*1];
        ham1.J_  = dim_ranges_[2*2];
        ham1.x_  = dim_ranges_[2*3];
        
        #pragma omp single
        {
          GlobalAttr = ham1.GetAttributes(); // assign attributes to GlobalAttr
          
          if (with_output)
            std::cout << "\n\nFPparams:\n" << FPparams << "\n"
                      << "*********************************************************" << "\n";
        }
    
        // Loop over values of the parameter space
        #pragma omp for collapse(4) schedule(dynamic,1)
        for (   int  t_idx=0;  t_idx<dim_lengths_[0]; ++t_idx)
         for (  int tp_idx=0; tp_idx<dim_lengths_[1]; ++tp_idx)
          for ( int  J_idx=0;  J_idx<dim_lengths_[2]; ++J_idx)
           for (int  x_idx=0;  x_idx<dim_lengths_[3]; ++x_idx)
           {
               // Adjust phase space parameters
               ham1.t_  = coord_arrays[0][t_idx];
               ham1.tp_ = coord_arrays[1][tp_idx];
               ham1.J_  = coord_arrays[2][J_idx];
               ham1.x_  = coord_arrays[3][x_idx]; // Assign value of x
               
               if (with_output) // Print current params
               {
                 std::cout << "\n\n";
                 if (dim_lengths_[0]!=1) std::cout <<  "t = " << ham1.t_  << "\t";
                 if (dim_lengths_[1]!=1) std::cout << "tp = " << ham1.tp_ << "\t";
                 if (dim_lengths_[2]!=1) std::cout <<  "J = " << ham1.J_  << "\t";
                 if (dim_lengths_[3]!=1) std::cout <<  "x = " << ham1.x_  << "\t";
                 std::cout << "\n";
               }
               
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
                   std::cout << "\tWARNING: failure to converge after limit reached"
                             <<  "\tt = " << ham1.t_
                             << "\ttp = " << ham1.tp_  
                             <<  "\tJ = " << ham1.J_ 
                             <<  "\tx = " << ham1.x_ << "\n";
                   ++numfails;
               }
               
               chi_s_grid  [idx(t_idx,tp_idx,J_idx,x_idx)] = ham1.MFs_.chi_s;
               chi_d_grid  [idx(t_idx,tp_idx,J_idx,x_idx)] = ham1.MFs_.chi_d;
               Delta_s_grid[idx(t_idx,tp_idx,J_idx,x_idx)] = ham1.MFs_.Delta_s;
               Delta_d_grid[idx(t_idx,tp_idx,J_idx,x_idx)] = ham1.MFs_.Delta_d;
               
               DeltaSC_s_grid   [idx(t_idx,tp_idx,J_idx,x_idx)] = DeltaSC_s;
               DeltaSC_d_grid   [idx(t_idx,tp_idx,J_idx,x_idx)] = DeltaSC_d;
               optweight_xx_grid[idx(t_idx,tp_idx,J_idx,x_idx)] = optweight_xx;
               optweight_yy_grid[idx(t_idx,tp_idx,J_idx,x_idx)] = optweight_yy;
               energy_grid      [idx(t_idx,tp_idx,J_idx,x_idx)] = energy;
               mu_grid          [idx(t_idx,tp_idx,J_idx,x_idx)] = mu;
               loops_grid       [idx(t_idx,tp_idx,J_idx,x_idx)] = loops; // Save the number of loops to array.
               
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


int pstudyA(const bool show_output)
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

int pstudyB(const bool show_output)
{
    // Order of coordinates: t, tp, J, x
    const size_t dim_lengths[4]  = {1, 101, 1, 31};
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

int pstudyB_fast(const bool show_output)
{
    // Order of coordinates: t, tp, J, x
    const size_t dim_lengths[4]  = {1, 26, 1, 16};
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


int main(int argc, char* argv[])
{
    const int info = pstudyA(false);
    return info;
}