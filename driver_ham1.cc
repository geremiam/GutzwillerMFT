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

// Classes that define parameter spaces for this Hamiltonian


class pspaceA_t { // Filling varied with interaction strength and temp held constant
  private:
    // Private copy constructor (prohibits copy creation)
    pspaceA_t(const pspaceA_t&);
    // Private assignment operator (prohibits assignment)
    const pspaceA_t& operator=(const pspaceA_t&);
    
  public:
    // x
    const size_t x_pts = 11;
    const double x_bounds [2] = {0.01, 0.31};
    
    const int parspace_pts = x_pts;
    
    // Coordinate variables
    double*const x_grid;// x coordinate variable
    
    // Variables for MF parameters. We store them in 1D arrays. The index keeps track of 
    // the point in parameter space (using the method index()).
    double*const chi_s_grid;
    double*const chi_d_grid;
    complex<double>*const Delta_s_grid;
    complex<double>*const Delta_d_grid;
    complex<double>*const DeltaSC_s_grid;
    complex<double>*const DeltaSC_d_grid;
    
    int   *const loops_grid; // holds the number of loops done at each point
    double*const energy_grid; // Holds the MF energy for later comparison
    double*const mu_grid; // Holds the chemical potential for the converged system
    
    // Constructor declaration
    pspaceA_t()
        :x_grid    (new double [x_pts]), // Coord vars
         chi_s_grid  (new double          [parspace_pts]), // Vars
         chi_d_grid  (new double          [parspace_pts]),
         Delta_s_grid(new complex<double> [parspace_pts]),
         Delta_d_grid(new complex<double> [parspace_pts]),
         DeltaSC_s_grid(new complex<double> [parspace_pts]),
         DeltaSC_d_grid(new complex<double> [parspace_pts]),
         loops_grid  (new int             [parspace_pts]),
         energy_grid (new double          [parspace_pts]),
         mu_grid     (new double          [parspace_pts]) // Important -- Note order
    {
        // The initialization list initializes the parameters and allocates memory for 
        // the arrays.
        
        // Coordinate variables are initialized
        const bool endpoint = true; // Include endpoint (not an important choice)
        LinInitArray(x_bounds[0], x_bounds[1], x_pts, x_grid, endpoint);
        
        // The other variables (the order parameters) span the coordinate space
        ValInitArray(parspace_pts, chi_s_grid,   -99.);//Init to unlikely value
        ValInitArray(parspace_pts, chi_d_grid,   -99.);//Init to unlikely value
        ValInitArray(parspace_pts, Delta_s_grid, -99.);//Init to unlikely value
        ValInitArray(parspace_pts, Delta_d_grid, -99.);//Init to unlikely value
        ValInitArray(parspace_pts, loops_grid,     -1); // Initialize to 0
        ValInitArray(parspace_pts, energy_grid,  -99.);
        ValInitArray(parspace_pts, mu_grid,       -9.);
        
        std::cout << "pspaceA_t instance created.\n";
    }
    // Destructor declaration
    ~pspaceA_t()
    {
        delete [] x_grid; // Coordinate variables
        delete [] chi_s_grid; // MF variables
        delete [] chi_d_grid;
        delete [] Delta_s_grid;
        delete [] Delta_d_grid;
        delete [] DeltaSC_s_grid;
        delete [] DeltaSC_d_grid;
        delete [] loops_grid;
        delete [] energy_grid;
        delete [] mu_grid;
        
        std::cout << "pspaceA_t instance deleted.\n";
    }
    
    int idx(const int f)
    {
        // Returns the index of the var arrays corresponding to indices i, j of the coord
        // var arrays. Trivial in this case.
        return f;
    }
    
    void SaveData(const string GlobalAttr, const string path)
    {
        // Method for saving the data of this class. This method uses the class from the 
        // module nc_IO that creates a simple NetCDF class and allows writing of variables.
        // We define parameters required to create the dataset. Don't forget to adjust these depending on the parameter space defined above. 
        // Important: Note that the order of the variables must be kept consistent.
        const size_t dims_num = 1;
        const string dim_names [dims_num] = {"x"};
        const size_t dim_lengths [dims_num] = {x_pts};
        
        const size_t vars_num = 6; // Variables other than coord variables
        string var_names [vars_num] = {"chi_s", "chi_d", "Delta_s", "Delta_d", "DeltaSC_s", "DeltaSC_d"}; // List for the variable names
        bool var_complex [vars_num] = {  false,   false,      true,      true,        true,        true}; // List for indicating whether vars are complex
        
        // Constructor for the dataset class creates a dataset
        newDS_t newDS(dims_num, dim_names, dim_lengths, vars_num, var_names, var_complex, GlobalAttr, path);
        
        newDS.DefCoordVar(0, "x"); // Define "default" coordinate variables
        
        newDS.EndDef(); // Exit definition mode
        
        newDS.WriteCoordVar(0, x_grid); // Write "default" coordinate variables
        
        // List for holding the pointers to the vars
        double* vars [vars_num] = {chi_s_grid, chi_d_grid, 
                                   reinterpret_cast<double*const>(Delta_s_grid), 
                                   reinterpret_cast<double*const>(Delta_d_grid),
                                   reinterpret_cast<double*const>(DeltaSC_s_grid), 
                                   reinterpret_cast<double*const>(DeltaSC_d_grid)};
        newDS.WriteVars(vars); // Write the variables
        
        newDS.WriteLoops(loops_grid); // Write loops variable
        newDS.WriteEnergy(energy_grid); // Write energy variable
        newDS.Writemu(mu_grid); // Write mu variable
    }
    
    int pstudy(const bool with_output = false)
    {
        // This routine performs the mean-field iterative search at every point in the 
        // parameter space defined above. Show output for diagnostics if with_output.
        int numfails = 0; // Tracks number of points which failed to converge after loops_lim
        
        MFs_t MFs_initial(0.1, 0., {0.,0.}, {0.1,0.});
    
        const double tol = 1.e-6;
        const int loops_lim = 800;
        const int mixing_vals_len = 7;
        const int   counter_vals [mixing_vals_len] = {100, 200, 300, 400, 500, 600, 700};
        const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3};
        FPparams_t FPparams(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    
        const int k1_pts = 298;//502; // Choose twice a prime number for grid resolution
        const int k2_pts = 298;//502;
        
        string GlobalAttr; // String in which to save the global attributes.
        
        // Must declare 'const' variables as firstprivate for compatibility with new gcc 
        // versions. See https://gcc.gnu.org/gcc-9/porting_to.html#ompdatasharing.
        #pragma omp parallel default(none) shared(FPparams,MFs_initial,GlobalAttr,std::cout) firstprivate(with_output) reduction(+:numfails)
        {
        // Declare and construct an instance of ham1_t
        ham1_t ham1(FPparams, MFs_initial, k1_pts, k2_pts);
    
        // Make any initial adjustment to the parameters
        ham1.set_nonzerotemp(5.e-3);
        ham1.t_  = 1.;
        ham1.tp_ = 0.;
        ham1.J_  = 1./5.;
        
        #pragma omp single
        {
          GlobalAttr = ham1.GetAttributes(); // assign attributes to GlobalAttr
          
          if (with_output)
            std::cout << "\n\nFPparams:\n" << FPparams << "\n"
                      << "*********************************************************" << "\n";
        }
    
        // Loop over values of the parameter space
        #pragma omp for schedule(dynamic,1)
        for (int f=0; f<x_pts; ++f)
        {
            if (with_output) // Print current params
              std::cout << "\n\nx = " << x_grid[f] << "\n";
            
            // Adjust phase space parameters
            ham1.x_ = x_grid[f]; // Assign value of x
            
            int loops=0; // Will be assigned the number of loops performed
            double mu=0.;// Will be assigned the chemical potential of the converged Ham.
            double energy=0.;
            complex<double> DeltaSC_s=0.;
            complex<double> DeltaSC_d=0.;
            const bool fail = ham1.FixedPoint(with_output, &loops, &mu, &energy, &DeltaSC_s, &DeltaSC_d);
            
            if (fail) // Print current params
            {
                std::cout << "\tWARNING: failure to converge after limit reached\t"
                          << "x = " << x_grid[f] << "\n";
                ++numfails;
            }
            
            chi_s_grid  [idx(f)] = ham1.MFs_.chi_s;
            chi_d_grid  [idx(f)] = ham1.MFs_.chi_d;
            Delta_s_grid[idx(f)] = ham1.MFs_.Delta_s;
            Delta_d_grid[idx(f)] = ham1.MFs_.Delta_d;
            
            DeltaSC_s_grid[idx(f)] = DeltaSC_s;
            DeltaSC_d_grid[idx(f)] = DeltaSC_d;
            energy_grid   [idx(f)] = energy;
            mu_grid       [idx(f)] = mu;
            loops_grid    [idx(f)] = loops; // Save the number of loops to array.
            
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

class pspaceB_t { // Filling varied with interaction strength and temp held constant
  private:
    // Private copy constructor (prohibits copy creation)
    pspaceB_t(const pspaceB_t&);
    // Private assignment operator (prohibits assignment)
    const pspaceB_t& operator=(const pspaceB_t&);
    
  public:
    // x
    const size_t x_pts = 31;
    const double x_bounds [2] = {0.01, 0.31};
    const size_t tp_pts = 101;
    const double tp_bounds[2] = {-0.5, 0.5};
    
    const int parspace_pts = x_pts*tp_pts;
    
    // Coordinate variables
    double*const x_grid; // x coordinate variable
    double*const tp_grid; // tp coordinate variable
    
    // Variables for MF parameters. We store them in 1D arrays. The index keeps track of 
    // the point in parameter space (using the method index()).
    double*const chi_s_grid;
    double*const chi_d_grid;
    complex<double>*const Delta_s_grid;
    complex<double>*const Delta_d_grid;
    complex<double>*const DeltaSC_s_grid;
    complex<double>*const DeltaSC_d_grid;
    
    int   *const loops_grid; // holds the number of loops done at each point
    double*const energy_grid; // Holds the MF energy for later comparison
    double*const mu_grid; // Holds the chemical potential for the converged system
    
    // Constructor declaration
    pspaceB_t()
        :x_grid    (new double [x_pts]), // Coord vars
         tp_grid   (new double [tp_pts]),
         chi_s_grid  (new double          [parspace_pts]), // Vars
         chi_d_grid  (new double          [parspace_pts]),
         Delta_s_grid(new complex<double> [parspace_pts]),
         Delta_d_grid(new complex<double> [parspace_pts]),
         DeltaSC_s_grid(new complex<double> [parspace_pts]),
         DeltaSC_d_grid(new complex<double> [parspace_pts]),
         loops_grid  (new int             [parspace_pts]),
         energy_grid (new double          [parspace_pts]),
         mu_grid     (new double          [parspace_pts]) // Important -- Note order
    {
        // The initialization list initializes the parameters and allocates memory for 
        // the arrays.
        
        // Coordinate variables are initialized
        const bool endpoint = true; // Include endpoint (not an important choice)
        LinInitArray(x_bounds[0], x_bounds[1], x_pts, x_grid, endpoint);
        LinInitArray(tp_bounds[0],tp_bounds[1],tp_pts,tp_grid,endpoint);
        
        // The other variables (the order parameters) span the coordinate space
        ValInitArray(parspace_pts, chi_s_grid,   -99.);//Init to unlikely value
        ValInitArray(parspace_pts, chi_d_grid,   -99.);//Init to unlikely value
        ValInitArray(parspace_pts, Delta_s_grid, -99.);//Init to unlikely value
        ValInitArray(parspace_pts, Delta_d_grid, -99.);//Init to unlikely value
        ValInitArray(parspace_pts, loops_grid,     -1); // Initialize to 0
        ValInitArray(parspace_pts, energy_grid,  -99.);
        ValInitArray(parspace_pts, mu_grid,       -9.);
        
        std::cout << "pspaceB_t instance created.\n";
    }
    // Destructor declaration
    ~pspaceB_t()
    {
        delete [] x_grid; // Coordinate variables
        delete [] tp_grid;
        delete [] chi_s_grid; // MF variables
        delete [] chi_d_grid;
        delete [] Delta_s_grid;
        delete [] Delta_d_grid;
        delete [] DeltaSC_s_grid;
        delete [] DeltaSC_d_grid;
        delete [] loops_grid;
        delete [] energy_grid;
        delete [] mu_grid;
        
        std::cout << "pspaceB_t instance deleted.\n";
    }
    
    int idx(const int tp_idx, const int x_idx)
    {
        // Returns the index of the var arrays corresponding to indices i, j of the coord
        // var arrays. Trivial in this case.
        return tp_idx*x_pts + x_idx;
    }
    
    void SaveData(const string GlobalAttr, const string path)
    {
        // Method for saving the data of this class. This method uses the class from the 
        // module nc_IO that creates a simple NetCDF class and allows writing of variables.
        // We define parameters required to create the dataset. Don't forget to adjust these depending on the parameter space defined above. 
        // Important: Note that the order of the variables must be kept consistent.
        const size_t dims_num = 2;
        const string dim_names [dims_num] = {"tp", "x"};
        const size_t dim_lengths [dims_num] = {tp_pts, x_pts};
        
        const size_t vars_num = 6; // Variables other than coord variables
        string var_names [vars_num] = {"chi_s", "chi_d", "Delta_s", "Delta_d", "DeltaSC_s", "DeltaSC_d"}; // List for the variable names
        bool var_complex [vars_num] = {  false,   false,      true,      true,        true,        true}; // List for indicating whether vars are complex
        
        // Constructor for the dataset class creates a dataset
        newDS_t newDS(dims_num, dim_names, dim_lengths, vars_num, var_names, var_complex, GlobalAttr, path);
        
        newDS.DefCoordVar(0, "tp");
        newDS.DefCoordVar(1, "x"); // Define "default" coordinate variables
        
        newDS.EndDef(); // Exit definition mode
        
        newDS.WriteCoordVar(0, tp_grid);
        newDS.WriteCoordVar(1, x_grid); // Write "default" coordinate variables
        
        // List for holding the pointers to the vars
        double* vars [vars_num] = {chi_s_grid, chi_d_grid, 
                                   reinterpret_cast<double*const>(Delta_s_grid), 
                                   reinterpret_cast<double*const>(Delta_d_grid),
                                   reinterpret_cast<double*const>(DeltaSC_s_grid), 
                                   reinterpret_cast<double*const>(DeltaSC_d_grid)};
        newDS.WriteVars(vars); // Write the variables
        
        newDS.WriteLoops(loops_grid); // Write loops variable
        newDS.WriteEnergy(energy_grid); // Write energy variable
        newDS.Writemu(mu_grid); // Write mu variable
    }
    
    int pstudy(const bool with_output = false)
    {
        // This routine performs the mean-field iterative search at every point in the 
        // parameter space defined above. Show output for diagnostics if with_output.
        int numfails = 0; // Tracks number of points which failed to converge after loops_lim
        
        MFs_t MFs_initial(0.1, 0., {0.,0.}, {0.1,0.});
    
        const double tol = 1.e-6;
        const int loops_lim = 800;
        const int mixing_vals_len = 7;
        const int   counter_vals [mixing_vals_len] = {100, 200, 300, 400, 500, 600, 700};
        const double mixing_vals [mixing_vals_len] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3};
        FPparams_t FPparams(tol, loops_lim, mixing_vals_len, counter_vals, mixing_vals);
    
        const int k1_pts = 298;//502; // Choose twice a prime number for grid resolution
        const int k2_pts = 298;//502;
        
        string GlobalAttr; // String in which to save the global attributes.
        
        // Must declare 'const' variables as firstprivate for compatibility with new gcc 
        // versions. See https://gcc.gnu.org/gcc-9/porting_to.html#ompdatasharing.
        #pragma omp parallel default(none) shared(FPparams,MFs_initial,GlobalAttr,std::cout) firstprivate(with_output) reduction(+:numfails)
        {
        // Declare and construct an instance of ham1_t
        ham1_t ham1(FPparams, MFs_initial, k1_pts, k2_pts);
    
        // Make any initial adjustment to the parameters
        ham1.set_nonzerotemp(5.e-3);
        ham1.t_  = 1.;
        ham1.J_  = 1./5.;
        
        #pragma omp single
        {
          GlobalAttr = ham1.GetAttributes(); // assign attributes to GlobalAttr
          
          if (with_output)
            std::cout << "\n\nFPparams:\n" << FPparams << "\n"
                      << "*********************************************************" << "\n";
        }
    
        // Loop over values of the parameter space
        #pragma omp for collapse(2) schedule(dynamic,1)
        for (int tp_idx=0; tp_idx<tp_pts; ++tp_idx)
         for (int x_idx=0;  x_idx<x_pts;  ++x_idx)
        {
            if (with_output) // Print current params
            {
              std::cout << "\n\ntp = " << tp_grid[tp_idx]
                        <<    "\tx = " <<  x_grid[x_idx]  << "\n";
            }
            
            // Adjust phase space parameters
            ham1.tp_ = tp_grid[tp_idx];
            ham1.x_  =  x_grid[x_idx]; // Assign value of x
            
            int loops=0; // Will be assigned the number of loops performed
            double mu=0.;// Will be assigned the chemical potential of the converged Ham.
            double energy=0.;
            complex<double> DeltaSC_s=0.;
            complex<double> DeltaSC_d=0.;
            const bool fail = ham1.FixedPoint(with_output, &loops, &mu, &energy, &DeltaSC_s, &DeltaSC_d);
            
            if (fail) // Print current params
            {
                std::cout << "\tWARNING: failure to converge after limit reached\t"
                          << "x = " << x_grid[x_idx] << "\ttp = " << tp_grid[tp_idx] << "\n";
                ++numfails;
            }
            
            chi_s_grid  [idx(tp_idx,x_idx)] = ham1.MFs_.chi_s;
            chi_d_grid  [idx(tp_idx,x_idx)] = ham1.MFs_.chi_d;
            Delta_s_grid[idx(tp_idx,x_idx)] = ham1.MFs_.Delta_s;
            Delta_d_grid[idx(tp_idx,x_idx)] = ham1.MFs_.Delta_d;
            
            DeltaSC_s_grid[idx(tp_idx,x_idx)] = DeltaSC_s;
            DeltaSC_d_grid[idx(tp_idx,x_idx)] = DeltaSC_d;
            energy_grid   [idx(tp_idx,x_idx)] = energy;
            mu_grid       [idx(tp_idx,x_idx)] = mu;
            loops_grid    [idx(tp_idx,x_idx)] = loops; // Save the number of loops to array.
            
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


int main(int argc, char* argv[])
{
    pspaceB_t pspace;
    int info = pspace.pstudy(false);
    
    return info;
}