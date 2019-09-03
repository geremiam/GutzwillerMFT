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
using std::complex;

class ham1_t
{
  /* Class for holding the Hamiltonian parameters. Parameters than are not likely to be 
  changed in a parameter study are made private (and const). Parameters that may be 
  varied during a parameter study are public, except for those on which other parameters 
  depend; such parameters are modifiable with a public method. */
  
  private:
    
    // Private copy constructor (prohibits copy creation)
    ham1_t(const ham1_t&);
    // Private assignment operator (prohibits assignment)
    const ham1_t& operator=(const ham1_t&);
    
  public:
    
    /* Parameters that the user doesn't need to modify after instantiation. */
    
    const double a_ = 1.; // a- and b-axis length
    const double c_ = 1.; // c-axis length
    
    int mixing_vals_len_;
    int*   counter_vals_;
    double* mixing_vals_;
    const int loops_lim_; // Limit to the number of iteration loops. Assigned in constructor.
    
    // Parameters for the momentum space grid. Assigned in the constructor.
    // Sets the number of momentum points IN THE FULL BZ
    const int ka_pts_; // Useful to allow the user to set these for convergence studies
    const int kb_pts_;
    const int kc_pts_;
    const int states_per_cell = 2; // Number of states in an (original) unit cell
    static const int num_harmonics = 4; // Number of harmonics considered.
    const int num_unit_cells = ka_pts_*kb_pts_*kc_pts_; // The number of (original) unit cells
    
    // We assume we are only working with harmonics given by G/2.
    const double Qa [num_harmonics] = {0., M_PI/a_, 0.,      M_PI/a_};
    const double Qb [num_harmonics] = {0., M_PI/a_, 0.,      M_PI/a_};
    const double Qc [num_harmonics] = {0., 0.,      M_PI/c_, M_PI/c_};
    
    const int addition_table [num_harmonics][num_harmonics]={{0, 1, 2, 3},
                                                             {1, 0, 3, 2},
                                                             {2, 3, 0, 1},
                                                             {3, 2, 1, 0}};
    
    // Size of the Hamiltonian, or equivalently number of bands
    const int num_bands = states_per_cell * num_harmonics; // The number of bands, i.e. the order of the matrix for each k
    const int ham_array_rows = num_bands; // Same as matrix order for full storage
    const int ham_array_cols = num_bands; // Same as matrix order for full storage
    const int num_states = num_unit_cells * states_per_cell;
    
    
    
    /* Parameters that have dependencies or are dependent on other parameters. These 
    should be modified with an appropriate method. */
    double rho_ = 0.4; // Average (global) electron density, between 0 and 1
    int filled_states=0; // Value assigned in constructor or upon calling 'assign_rho()'
    
    bool zerotemp_ = true; // Sets whether the temperature is zero or not
    double T_ = 0.; // Sets the temperature. Irrelevant when zerotemp_ is true.
    
    // Function for performing the direct product in the indices Q and alpha.
    int idx(const int Q, const int alpha) const;
    
    
    // These functions are used in the method ComputeMFs()
    void AddContribution_rho(const double*const occs, const complex<double>*const*const evecs, double*const rho_A, double*const rho_B) const;
    
    void Assign_h(double ka, double kb, double kc, complex<double>*const h) const;
    void Assign_V(const int Q, complex<double>*const V) const;
    //void Assign_V_manual(const int Q, complex<double>*const V) const;
    void Assign_ham(const double ka, const double kb, const double kc, complex<double>*const*const ham_array) const;
    
    void ComputeMFs_old(double*const rho_s_out, double*const rho_a_out, double*const HFE_p=NULL, double*const mu_p=NULL) const;
    void ComputeMFs    (double*const rho_s_out, double*const rho_a_out, double*const HFE_p=NULL, double*const mu_p=NULL) const;
    
    
  //public:
    
    /* Settings for the iterative search */
    
    // The MF values that are used in the Hamiltonian.
    // This is where the user sets the initial values of the MFs.
    
    // Because each harmonic is given by some G/2, the densities are all real.
    double rho_s_ [num_harmonics] = {0.};
    double rho_a_ [num_harmonics] = {0.};
    
    double HFE_ = -99.; // For storing the free energy
    double mu_  = -9.;  // For storing the chemical potential
    
    const double tol_; // Tolerance for the equality of the mean fields
    
    
    // Hamiltonian parameters that the user may want to change
    double t1_  = 1.; // x- and y-direction hopping
    double t1p_ = 0.7; // a- and b-direction hopping
    double t2A_ = 0.1; // "A-emanating" solenoid-like hopping
    double t2B_ = 0.4; // "B-emanating" solenoid-like hopping
    double t3_  = 0.5; // c-direction hopping
    double V1_  = 1.0; // Repulsion between neighbours in x and y directions
    double V1p_ = 0.7; // Repulsion between neighbours in a and b directions
    double V2_  = 0.4; // We choose to use the same V2 for A- and B-emanating
    double V3_  = 0.5; // Repulsion between neighbours in c direction
    
    
    void resetMFs(); // Resets MFs to default starting values.
    
    void assign_rho(const double rho); // Assign rho_ and dependent variables
    void set_zerotemp();
    void set_nonzerotemp(const double T);
    
    // Constructor declaration
    ham1_t(const int ka_pts=62, const int kb_pts=62, const int kc_pts=62, const double tol=4.e-6, const int loops_lim=300,
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
