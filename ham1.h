// ham1.h
/* Description.
   Here is what the user needs to change upon changing the harmonics: 
   (1) the arrays Qa, Qb, and Qc;
   (2) the array addition_table; 
   (3) the implementation of Assign_V_manual. */
#ifndef HAM1_H
#define HAM1_H

#include <complex>
#include <cmath> // For the constant M_PI
#include <ostream>
#include <string>
#include "kspace.h" // Defines a class for holding a band structure
using std::complex;
using std::string;

class MFs_t
{
// Since this class contains only Plain Old Data types, it is not necessary to program a
// copy constructor or a copy assignment operator.
  private:
    string output_text_; // Used in const char* casting
  public:
    complex<double> chi_s = {0.,0.};
    complex<double> chi_d = {0.,0.};
    complex<double> Delta_s = {0.,0.};
    complex<double> Delta_d = {0.,0.};
    
    
    // Constructor declarations
    MFs_t(complex<double> chi_s_=0., complex<double> chi_d_=0., complex<double> Delta_s_=0., complex<double> Delta_d_=0.);
    MFs_t(const complex<double>*const MFs_array);
    
    void MFs_from_array(const complex<double>*const MFs_array);
    void MFs_to_array(complex<double>*const MFs_array) const;
    
    void update_values(const MFs_t& MFs_new, const double mixratio);
    
    bool check_bound(const double bound) const;
    
    // Operator for casting to const char*.
    operator const char* ();
    
    // Binary addition operator
    MFs_t operator + (const MFs_t& MFs_toadd);
    
    // Multiplication by a double
    MFs_t operator * (const double factor);
    
    // Negation operator
    MFs_t operator - ();
};

class ham1_t
{
  /* Class for holding the Hamiltonian parameters. Parameters that are not likely to be 
  changed in a parameter study are made private (and const). Parameters that may be 
  varied during a parameter study are public, except for those on which other parameters 
  depend; such parameters are modifiable with a public method. */
  
  private:
    
    // Private copy constructor (prohibits copy creation)
    ham1_t(const ham1_t&);
    // Private assignment operator (prohibits assignment)
    const ham1_t& operator=(const ham1_t&);
    
    
    /* Parameters that the user doesn't need to modify after instantiation. */
    
    int mixing_vals_len_;
    int*   counter_vals_;
    double* mixing_vals_;
    const int loops_lim_; // Limit to the number of iteration loops. Assigned in constructor.
    
    // Parameters for the momentum space grid. Assigned in the constructor.
    const double b1_[2] = {1., 0.};
    const double b2_[2] = {0., 1.};
    const int k1_pts_; // Useful to allow the user to set these for convergence studies
    const int k2_pts_;
    const int states_per_cell = 2; // Number of states in an (original) unit cell
    const int num_unit_cells = k1_pts_*k2_pts_; // The number of (original) unit cells
    
    // Size of the Hamiltonian, or equivalently number of bands
    const int num_bands = states_per_cell; // The number of bands, i.e. the order of the matrix for each k
    const int ham_array_rows = num_bands; // Same as matrix order for full storage
    const int ham_array_cols = num_bands; // Same as matrix order for full storage
    const int num_states = num_unit_cells * states_per_cell;
    
    double bisec1(const double a_in, const double b_in, const bool show_output=false) const;
    void diag(const double kx, const double ky, const double mu_local, double& E_cal, double& u, complex<double>& v) const;
    double utility(const double mu_local) const;
    void chempot();
    MFs_t compute_MFs(const double kx, const double ky) const;
    
    
    //void ComputeMFs_old(double*const rho_s_out, double*const rho_a_out, double*const HFE_p=NULL, double*const mu_p=NULL) const;
    //void ComputeMFs    (double*const rho_s_out, double*const rho_a_out, double*const HFE_p=NULL, double*const mu_p=NULL) const;
    
    kspace_t kspace; // Initialized in the initializer list.
    
    double HFE_ = -99.; // For storing the free energy
    double mu_  = -9.;  // For storing the chemical potential
    
    const double tol_; // Tolerance for the equality of the mean fields
    
  public:
    
    /* Settings for the iterative search */
    
    // The MF values that are used in the Hamiltonian.
    // This is where the user sets the initial values of the MFs.
    
    // Declare instances of the MFs_t class, which holds the MFs.
    MFs_t MFs_initial;
    MFs_t MFs;
    
    
    // Hamiltonian parameters that the user may want to change
    double T_ = 1.; // Sets the temperature.
    double x_ = 0.1; // Hole doping; between -1 and 1.
    double t_  = 1.; // NN hopping amplitude
    double tp_ = -0.25; // NNN hopping amplitude
    double J_  = 1./3.; // Heisenberg repulsion
    
    // Resets MFs to default starting values.
    void resetMFs();
    
    
    // Constructor declaration
    ham1_t(const int k1_pts=62, const int k2_pts=62, const double tol=4.e-6, const int loops_lim=300,
           const int mixing_vals_len=0, const int*const counter_vals=NULL, const double*const mixing_vals=NULL);
    ~ham1_t(); // Destructor declaration
    
    bool FixedPoint(int*const num_loops_p=NULL, const bool with_output=false);
    
    std::string GetAttributes();
    
    // ROUTINES FOR CALCULATING THE FREE ENERGY
    double Helmholtz  (const double*const energies, const double mu, const double*const rho_s_out, const double*const rho_a_out) const;
    double Omega_trial(const double*const energies, const double mu, const double*const rho_s_out, const double*const rho_a_out) const;
    double Omega_MF   (const double*const energies, const double mu) const;
    double mean_V   (const double*const rho_s_out, const double*const rho_a_out) const;
    double mean_V_MF(const double*const rho_s_out, const double*const rho_a_out) const;
};

#endif
