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
    double            chi_s = 0.;
    double            chi_d = 0.;
    complex<double> Delta_s = {0.,0.};
    complex<double> Delta_d = {0.,0.};
    
    
    // Constructor declarations
    MFs_t(double chi_s_=0., double chi_d_=0., complex<double> Delta_s_=0., complex<double> Delta_d_=0.);
    
    bool check_bound(const double bound) const;
    
    // Gives output format
    string output_format() const;
    
    // Operator for casting to const char*.
    operator const char* ();
    
    // Binary addition operator
    MFs_t operator + (const MFs_t& MFs_toadd);
    
    // Binary subtraction operator
    MFs_t operator - (const MFs_t& MFs_toadd);
    
    // Multiplication by a double
    MFs_t operator * (const double factor);
    
    // Negation operator
    MFs_t operator - ();
};

class FPparams_t
{
  private:
    const FPparams_t& operator = (const FPparams_t&); // Private assignment operator
    
    string output_text_; // Used in const char* casting
    
    const int mixing_vals_len_;
    int*         counter_vals_;
    double*       mixing_vals_;
    
  public:
    const double    tol_;
    const int loops_lim_;
    FPparams_t(const double tol, const int loops_lim, const int mixing_vals_len, const int*const counter_vals, const double*const mixing_vals);
    ~FPparams_t();
    FPparams_t(const FPparams_t& copy_source); // Copy constructor
    
    double mixratio(const int iteration, const bool silent=false) const;
    operator const char* (); // Operator for casting to const char*
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
    
    
    // Parameters that the user doesn't need to modify after instantiation.
    
    // Parameters for the momentum space grid. Assigned in the constructor.
    const double b1_[2] = {2.*M_PI, 0.     };
    const double b2_[2] = {0.,      2.*M_PI};
    const int k1_pts_; // Useful to allow the user to set these for convergence studies
    const int k2_pts_;
    const int states_per_cell = 2; // Number of states in an (original) unit cell
    const int num_unit_cells = k1_pts_*k2_pts_; // The number of (original) unit cells
    const int num_states = num_unit_cells * states_per_cell;
    const int num_bands = states_per_cell; // Size of the Hamiltonian, or equivalently number of bands
    
    // Declare an instace of MFs_t to hold the starting MF values.
    MFs_t MFs_initial_;
    
    kspace_t kspace; // Initialized in the initializer list.
    
    FPparams_t FPparams_; // Parameters for the fixed-point algorithm. Assigned in constructor.
    
    // Parameters that the user must modify via methods because they are interdependent
    bool zerotemp_ = false;
    double T_ = 0.1; // Temperature
    
    // Useful internal methods
    double disp(const double kx, const double ky) const;
    double g_t() const;
    double g_S() const;
    
    bool   diag(const double kx, const double ky, const double mu_local, double& E, double& u, complex<double>& v) const;
    double chempot_utility(const double mu_local) const;
    double bisec1(const double a_in, const double b_in, const bool show_output=false) const;
    double chempot(const bool show_output=false) const;
    MFs_t  compute_MFs(double*const mu_output=NULL, double*const kin_energy_p=NULL, double*const optweight_xx_p=NULL, double*const optweight_yy_p=NULL) const;
        
  public:
    // Constructor and destructor declarations
    ham1_t(const FPparams_t FPparams, const MFs_t MFs_initial, const int k1_pts, const int k2_pts);
    ~ham1_t();
    
    // Hamiltonian parameters that the user may want to change
    double x_ = 0.1; // Hole doping; between -1 and 1.
    double t_  = 1.; // NN hopping amplitude
    double tp_ = -0.25; // NNN hopping amplitude
    double J_  = 1./3.; // Heisenberg repulsion
    // Methods for modifying interdependent parameters
    void set_zerotemp();
    void set_nonzerotemp(const double T);
    
    MFs_t MFs_; // Used to store the current MF values
    
    bool FixedPoint(const bool with_output=false, int*const num_loops_p=NULL, double*const mu_output=NULL, double*const energy_p=NULL, 
                    complex<double>*const DeltaSC_s=NULL, complex<double>*const DeltaSC_d=NULL,
                    double*const optweight_xx_p=NULL, double*const optweight_yy_p=NULL);
    std::string GetAttributes();
};

#endif
