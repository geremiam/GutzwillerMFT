// ham1.cc
/* Description */

#include <iostream>
#include <iomanip> // For the function std::setprecision()
#include <sstream> // For stringstreams
#include <complex> // Needed to import alloc.h
#include <cassert> // For assert function
#include <cmath>
#include "ticktock.h"
#include "alloc.h" // Allocation/deallocation of arrays
#include "array_init.h" // Initialization of arrays
#include "math.h" // Various math functions
//#include "misc_routines.h"
//#include "chempot.h"
#include "kspace.h" // Defines a class for holding a band structure
#include "diag.h" // Routines for finding evals and evecs
#include "ham1.h" // Include header file for consistency check
//using std::cos; using std::sin; using std::conj; // Not sure if these are necessary
using std::polar;
using std::to_string;
using std::complex;
using std::abs;


// **************************************************************************************
// Implementation of the class MFs_t
MFs_t::MFs_t(double chi_s_, double chi_d_, complex<double> Delta_s_, complex<double> Delta_d_)
    :chi_s(chi_s_), chi_d(chi_d_), Delta_s(Delta_s_), Delta_d(Delta_d_)
{}
bool MFs_t::check_bound(const double bound) const
{
    if (bound<=0.)
        std::cout << "\n*** WARNING: argument of method 'MFs_t::check_bound()' should be greater than zero.\n";
    return (abs(chi_s)<bound) && (abs(chi_d)<bound) && (abs(Delta_s)<bound) && (abs(Delta_d)<bound);
}
string MFs_t::output_format() const
{
    std::ostringstream output;
    output << "chi_s    chi_d    Delta_s    Delta_d";
    return output.str();
}
MFs_t::operator const char* () // Allows output as a string stream.
{
    std::ostringstream output;
    output << std::scientific << std::showpos << chi_s << "\t" << chi_d << "\t" << Delta_s << "\t" << Delta_d;
    output_text_ = output.str();
    return output_text_.c_str();
}
MFs_t MFs_t::operator + (const MFs_t& MFs_toadd) // Addition of two instances
{
    MFs_t newMFs(chi_s+MFs_toadd.chi_s, chi_d+MFs_toadd.chi_d, Delta_s+MFs_toadd.Delta_s, Delta_d+MFs_toadd.Delta_d);
    return newMFs;
}
MFs_t MFs_t::operator - (const MFs_t& MFs_toadd) // Subtraction of one instance from another
{
    //MFs_t newMFs = *this + (-MFs_toadd);
    MFs_t newMFs(chi_s-MFs_toadd.chi_s, chi_d-MFs_toadd.chi_d, Delta_s-MFs_toadd.Delta_s, Delta_d-MFs_toadd.Delta_d);
    return newMFs;
}
MFs_t MFs_t::operator * (const double factor) // Multiplication by a double
{
    MFs_t newMFs(factor*chi_s, factor*chi_d, factor*Delta_s, factor*Delta_d);
    return newMFs;
}
MFs_t MFs_t::operator - () // Negation
{
    MFs_t newMFs = *this * -1.; // Invokes copy assignment operator
    return newMFs;
}

// **************************************************************************************
// Implementation of the class FPparams_t
FPparams_t::FPparams_t(const double tol, const int loops_lim, const int mixing_vals_len, const int*const counter_vals, const double*const mixing_vals)
    :tol_(tol), loops_lim_(loops_lim), mixing_vals_len_(mixing_vals_len)
{
    // Check some basic stuff.
    assert(counter_vals!=NULL);
    assert(mixing_vals!=NULL);
    assert(mixing_vals_len>0);
    assert(loops_lim>0);
    assert(tol>0.);
    
    counter_vals_ = new int    [mixing_vals_len_];
    mixing_vals_  = new double [mixing_vals_len_];
    
    for (int i=0; i<mixing_vals_len_; ++i)
    {
      counter_vals_[i] = counter_vals[i];
      mixing_vals_[i]  = mixing_vals[i];
    }
    
    // Check that values provided are strictly increasing.
    for (int i=0; i<mixing_vals_len_-1; ++i)
    {
      assert(counter_vals_[i]<counter_vals_[i+1]);
      assert(mixing_vals_[i]>mixing_vals_[i+1]);
    }
}
FPparams_t::~FPparams_t()
{
    delete [] counter_vals_;
    delete [] mixing_vals_;
}
FPparams_t::FPparams_t(const FPparams_t& copy_source) // Copy constructor
    :tol_(copy_source.tol_), loops_lim_(copy_source.loops_lim_), mixing_vals_len_(copy_source.mixing_vals_len_)
{
    counter_vals_ = new int    [mixing_vals_len_];
    mixing_vals_  = new double [mixing_vals_len_];
    
    for (int i=0; i<mixing_vals_len_; ++i)
    {
      counter_vals_[i] = copy_source.counter_vals_[i];
      mixing_vals_[i]  = copy_source.mixing_vals_[i];
    }
}
double FPparams_t::mixratio(const int counter, const bool silent) const
{
    /* When 'counter' reaches one of the values in 'counter_values', the corresponding 
    entry of 'mixratio_values' is returned. Values in counter_vals must be strictly increasing.*/
    double mixratio = 1.; // Starting value should be 1.
    for (int i=0; i<mixing_vals_len_; ++i)
    {
        if (counter>=counter_vals_[i]) 
            mixratio = mixing_vals_[i];
        if (counter==counter_vals_[i]) 
        {
            if (!silent)
                std::cout << "\n\t Counter has reached " << counter
                          << ".\tmixratio = " << mixratio << "\n\n";
        }
    }
    
    return mixratio;
}
FPparams_t::operator const char* ()
{
    std::ostringstream output;
    output << "tol_ = " << tol_ << ", loops_lim_ = " << loops_lim_ << ", mixing_vals_len_ = " << mixing_vals_len_ << ",";
    output << "\ncounter_vals_ = \n";
    for (int i=0; i<mixing_vals_len_; ++i)
        output << counter_vals_[i] << " ";
    output << "\nmixing_vals_ = \n";
    for (int i=0; i<mixing_vals_len_; ++i)
        output << mixing_vals_[i] << " ";
    output_text_ = output.str();
    return output_text_.c_str();
}
// **************************************************************************************
//Implementation of the class ham1_t

double ham1_t::chempot_utility(const double mu_local) const
{
    // Function whose zero gives the "correct" value of mu.
    // Remember, std::norm() is square magnitude.
    double ans = 0.;
    
    // Tracks the number of marginal cases within the BZ, where we set |u|^2 = |v|^2 = 1/2.
    int marginals = 0;
    
    //#pragma omp parallel default(none) firstprivate(mu_local) shared(T_,num_unit_cells,k1_pts_,k2_pts_,kspace) reduction(+:ans, marginals)
    {
    double          E=0.;
    double          u=0.;
    complex<double> v=0.;
    //#pragma omp for collapse(2)
    for (int i=0; i<k1_pts_; ++i)
      for (int j=0; j<k2_pts_; ++j)
      {
        const double kx = kspace.kx_grid[kspace.k_i(i,j)];
        const double ky = kspace.ky_grid[kspace.k_i(i,j)];
        marginals += diag(kx, ky, mu_local, E, u, v); // Assign values to u, v, and E.
        
        double factor=0.;
        if (zerotemp_)
            factor = 1.;
        else
            factor = 1. - 2.*nF(T_,E);
        
        ans += (std::norm(u) - std::norm(v)) * factor / num_unit_cells;
      }
    }
    ans -= x_;
    
    if (marginals!=0)
        std::cout << "WARNING: number of marginal cases is " << marginals << " out of " 
        << num_unit_cells << ", i.e. " << 100.*marginals/num_unit_cells << "%." << std::endl;
    
    return ans;
}

double ham1_t::bisec1(const double a_in, const double b_in, const bool show_output) const
{
    // The chemical potential is found using the bisection method. Uses 100 times machine 
    // epsilon as relative tolerance.
    
    // Make sure that a < b.
    if ( !(a_in<b_in) )
      std::cout << "WARNING: function \"bisec\" requires a_in<b_in." << std::endl;
    // Make sure there is a zero between the two inputs.
    if ( !(chempot_utility(a_in)*chempot_utility(b_in)<0.) )
    {
      std::cout << "\nWARNING: the function does not change sign between a_in and b_in." << std::endl;
      std::cout << "a = " << a_in << ",\tb = " << b_in << "\n";
      std::cout << "f(a) = " << chempot_utility(a_in) << ",\tf(b) = " << chempot_utility(b_in) << "\n";
    }
    const bool increasing = (chempot_utility(a_in) < chempot_utility(b_in)); // Determine if the function is increasing
    
    // Starting values for a and b. Choose them slightly outside the energy range to be safe.
    double a = a_in;
    double b = b_in;
    
    int counter = 0;
    bool converged = false;
    
    do // Update via Newton's method
    {
      counter++; // Increment 'counter' by one
      
      const double midpoint = (a+b)/2.; // Get the midpoint between a and b
      // Find the image of the function at the midpoint
      const double image = chempot_utility(midpoint);
      
      if (show_output)
          std::cout << "a = " << a << ", b = " << b << ", b-a = " << b-a << "\t"
                    << "midpoint = " << midpoint << ", image = " << image;
      
      if (image==0.)
      {
        a = midpoint; // If image is 0, the loop terminates and the value is returned.
        b = midpoint;
      }
        
      else if (image<0.)
      {
        if (increasing)
          a = midpoint; // If 'image' is negative, 'midpoint' is the new 'a'
        else
          b = midpoint;
      }
      else if (image>0.)
      {
        if (increasing)
          b = midpoint; // If 'image' is positive, 'midpoint' is the new 'b'
        else
          a = midpoint;
      }
      else
        std::cout << "\nWARNING: image has invalid value." << std::endl;
      
      const double tolerance = 100. * eps_double * std::abs((a+b))/2.;//eps gives rel. tol.
      if (show_output)
          std::cout << "\ttolerance = " << tolerance << std::endl;
      if (std::abs((b-a))<tolerance)
        converged = true; // Check if the half difference is below tolerance
      
    } while (!converged); // The looping stops if convergence has been reached.
    
    return (b+a)/2.; // Take mu to be the average of a and b.
}

double ham1_t::chempot(const bool show_output) const
{
    // Finds the chemical potential using a bisection method. To be called before computing the MFs.
    const double max = 5. * (abs(t_) + abs(tp_) + abs(J_));
    const double mu_out = bisec1(-max, max, show_output);
    return mu_out;
}

bool ham1_t::diag(const double kx, const double ky, const double mu_local, double& E, double& u, complex<double>& v) const
{
    // Calculates three intermediate quantities at a given momentum for the current 
    // parameter values. The correct chemical potential must be given as an argument.
    
    bool marginal = false; // Signals marginal cases, where we set |u|^2 = |v|^2 = 1/2.
    
    const complex<double>& chi_s_ = MFs_.chi_s;
    const complex<double>& chi_d_ = MFs_.chi_d;
    const complex<double>& Delta_s_ = MFs_.Delta_s;
    const complex<double>& Delta_d_ = MFs_.Delta_d;
    
    // Expressions appearing in the Bloch Hamiltonian
    // ** Note that, in our case, eps_p and eps_m are one and the same. **
    const double        eps = -2.*t_* (cos(kx)+cos(ky)) - 4.*tp_*cos(kx)*cos(ky);
    const double          f = -3./2.*J_*std::real( cos(kx)*(chi_s_+chi_d_)       + cos(ky)*(chi_s_-chi_d_) );
    const complex<double> g =  3./2.*J_*         ( cos(kx) * (Delta_s_+Delta_d_) + cos(ky) * (Delta_s_-Delta_d_) );
    
    // Renormalization factors
    const double g_t = 2.*x_/(1.+x_);
    const double g_S = 4./((1.+x_)*(1.+x_));
    
    const double         xi = g_t*eps + g_S*f - mu_local;
    const complex<double> D = g_S*g;
    
    E = sqrt( xi*xi + std::norm(D) ); // Remember, std::norm() gives magnitude squared.
    
    // Carefully handle 0./0. situations.
    if ( (xi==0.) && (E==0.) )
    {
        u =                      std::sqrt(1./2.);
        v = - polar(1.,arg(D)) * std::sqrt(1./2.); // Expensive to find the phase factor.
        marginal = true;
    }
    else
    {
        u =                      std::sqrt( (1.+xi/E)/2. );
        v = - polar(1.,arg(D)) * std::sqrt( (1.-xi/E)/2. );
    }
    
    // Check for nans.
    if (std::isnan(E))
        std::cout << "WARNING: E is nan. mu_local = " << mu_local << std::endl;
    if (std::isnan(u))
        std::cout << "WARNING: u is nan. mu_local = " << mu_local << std::endl;
    if ( std::isnan(v.real()) || std::isnan(v.imag()) )
        std::cout << "WARNING: v is nan. mu_local = " << mu_local << std::endl;
    
    return marginal;
}

MFs_t ham1_t::compute_MFs(double*const mu_output_p) const
{
    // Step 1: calculate chemical potential for these parameters
    const double mu_local = chempot();
    if (mu_output_p!=NULL)
        *mu_output_p = mu_local; // Assign the computed value to mu_output_p
    
    // Step 2: calculate MFs using chemical potential
    // To be called after the correct mu is found and assigned.
    
    // Tracks the number of marginal cases within the BZ, where we set |u|^2 = |v|^2 = 1/2.
    int marginals = 0;
    
    // Variables for reduction
    double                x = 0.;
    double            chi_x = 0.;
    double            chi_y = 0.;
    complex<double> Delta_x = 0.;
    complex<double> Delta_y = 0.;
    
    //#pragma omp declare reduction(+:complex<double>:omp_out+=omp_in) // Must declare reduction on complex numbers
    //#pragma omp parallel default(none) firstprivate(mu_local) shared(T_,num_unit_cells,k1_pts_,k2_pts_,kspace) reduction(+:x, chi_x, chi_y, Delta_x, Delta_y, marginals)
    {
    //#pragma omp for collapse(2)
    for (int i=0; i<k1_pts_; ++i)
      for (int j=0; j<k2_pts_; ++j)
      {
        // Variables to be assigned to
        double          E=0.;
        double          u=0.;
        complex<double> v=0.;
        
        const double kx = kspace.kx_grid[kspace.k_i(i,j)];
        const double ky = kspace.ky_grid[kspace.k_i(i,j)];
        marginals += diag(kx, ky, mu_local, E, u, v); // Assign values to u, v, and E.
        
        double factor=0.;
        if (zerotemp_)
            factor = 1. / num_unit_cells;
        else
            factor = (1.-2.*nF(T_, E)) / num_unit_cells;
        
        x       +=         factor *           (std::norm(u)-std::norm(v));
        chi_x   += - 0.5 * factor * cos(kx) * (std::norm(u)-std::norm(v));
        chi_y   += - 0.5 * factor * cos(ky) * (std::norm(u)-std::norm(v));
        Delta_x += -       factor * cos(kx) * u * v; // Needs to be tested
        Delta_y += -       factor * cos(ky) * u * v; // Needs to be tested
      }
    }
    
    if (abs(x-x_)>1.e-14)
        std::cout << "\nWARNING: output holedoping does not match input holedoping."
                  << "\nx-x_ = " << x-x_ << std::endl;
    
    if (marginals!=0)
        std::cout << "WARNING: number of marginal cases is " << marginals << " out of " 
        << num_unit_cells << ", i.e. " << 100.*marginals/num_unit_cells << "%." << std::endl;
    
    const double            chi_s = (chi_x  +  chi_y) / 2.;
    const double            chi_d = (chi_x  -  chi_y) / 2.;
    const complex<double> Delta_s = (Delta_x + Delta_y)/2.;
    const complex<double> Delta_d = (Delta_x - Delta_y)/2.;
    
    MFs_t MFs_out(chi_s, chi_d, Delta_s, Delta_d);
    
    return MFs_out;
}

bool ham1_t::FixedPoint(const bool with_output, int*const num_loops_p, double*const mu_output_p, double*const energy_p, complex<double>*const DeltaSC_s, complex<double>*const DeltaSC_d)
{
    // Performs the iterative self-consistent search using the parameters from ham1. 
    // The initial values of the MF attributes are used as the starting values for the 
    // search; the end values are also stored in these same MF attributes.
    if (with_output) // Format display output and set precision
        std::cout << std::scientific << std::showpos;
    
    // Declare MFs variables and INITIALIZE THEM TO INPUT VALUES
    MFs_ = MFs_initial_;
    MFs_t MFs_out(MFs_initial_);
    double energy=0., energy_prev=0.; // For keeping track of free energy
    
    if (with_output)
       std::cout << MFs_out.output_format() << std::endl; // Order of MFs in output
    
    int counter = 0; // Define counter for number of loops
    bool converged=false, fail=false; // Used to stop the while looping
    do // Iterate until self-consistency is achieved
    {
        ++counter; // Increment counter
        
        // Past a certain number of loops, we mix in part of the previous input vals
        // Mixing fraction mixratio (mixratio=1 corresponds to using fully new value)
        const double mixratio = FPparams_.mixratio(counter, !with_output);
        
        MFs_ = MFs_*(1.-mixratio) + MFs_out*mixratio; // Update mean-field values
        energy_prev = energy; // Store previous free energy in energy_prev
        
        
        if (with_output)
        {
          std::cout << counter << "  mixratio = " << mixratio << "\ninput\t";
          std::cout << MFs_;
        }
        
        // Compute MFs. Output is assigned to MFs_out.
        // The chemical potential is stored in *mu_output_p. In the end, the chemical potential for the converged MFs is output.
        //ComputeMFs(rho_s_out, rho_a_out, &HFE_, &mu_);
        MFs_out = compute_MFs(mu_output_p);
        
        if (with_output)
          std::cout << "|  " << energy << std::endl;
        
        MFs_t MFs_diff = MFs_out - MFs_; // Differences between outputs and inputs
        
        if (with_output) // Print the differences
        {
          std::cout << "diff\t";
          std::cout << MFs_diff;
          std::cout << "|  " << energy-energy_prev << std::endl;
        }
        
        // Test for convergence.
        converged = MFs_diff.check_bound(FPparams_.tol_);
        fail = (!converged) && (counter>FPparams_.loops_lim_); // Must come after converged line
        
    } while (!converged && !fail);
    
    // Unless num_loops_p is the null pointer, assign the number of loops to its location
    if (num_loops_p!=NULL) *num_loops_p = counter;
    // The free energy for the converged MFs is assigned to *energy_p.
    if (energy_p!=NULL) *energy_p = energy;
    // Calculate the full SC order parameter and assign it, if applicable
    const double g_t = 2.*x_/(1.+x_);
    if (DeltaSC_s!=NULL) *DeltaSC_s = g_t * MFs_.Delta_s;
    if (DeltaSC_d!=NULL) *DeltaSC_d = g_t * MFs_.Delta_d;
    
    // We make sure that one of "converged" or "fail" is true.
    assert(converged != fail);
    
    return fail;
}




// Constructor and destructor
ham1_t::ham1_t(const FPparams_t FPparams, const MFs_t MFs_initial, const int k1_pts, const int k2_pts)
    :FPparams_(FPparams), MFs_initial_(MFs_initial),
    k1_pts_(k1_pts), k2_pts_(k2_pts), 
    kspace(b1_, b2_, k1_pts, k2_pts) // We do not need memory for evals or evecs
{
    std::cout << "ham1_t instance created.\n";
}

ham1_t::~ham1_t()
{
    /* Destructor implementation */
    std::cout << "ham1_t instance deleted.\n";
}


std::string ham1_t::GetAttributes()
{
    // Define a string of metadata
    std::ostringstream strs; // Declare a stringstream to which to write the attributes
    strs << "**ham1 parameters**"
         << "\nk-space points: " << k1_pts_ << "," << k1_pts_ << "; num_bands = " << num_bands
         << ", zerotemp = " << zerotemp_ << ", T = " << T_ << "\nx = " << x_
         << "\nt = " << t_ << ", tp = " << tp_ << ", J = " << J_
         << "\n" << FPparams_
         << "\nStarting values of the MFs: " << MFs_initial_;
        
    const std::string Attributes = strs.str(); // Write the stream contents to a string
    
    return Attributes;
}

void ham1_t::set_zerotemp()
{
    zerotemp_ = true;
    T_ = -1.;
}
void ham1_t::set_nonzerotemp(const double T)
{
    zerotemp_ = false;
    T_ = T;
    if (T<=0.)
        std::cout << "WARNING: set_nonzerotemp() requires a strictly positive argument.\n";
}

// 
// // **************************************************************************************
// // ROUTINES FOR CALCULATING THE FREE ENERGY
// 
// double ham1_t::Helmholtz(const double*const energies, const double mu, const double*const rho_s_out, const double*const rho_a_out) const
// {
//     // The Helmholtz FE normalized according to the NUMBER OF (ORIGINAL) UNIT CELLS
//     // Should be used to compare states of a fixed particle number
//     // Second term is + mu n (properly normalized).
//     // We use the "real" density rho_s_[0] instead of the target density rho_. These may 
//     // be different because of filled_states gets rounded to an int.
//     const double ans = Omega_trial(energies,mu,rho_s_out,rho_a_out) + mu*rho_s_[0]*states_per_cell;
//     return ans;
// }
// 
// double ham1_t::Omega_trial(const double*const energies, const double mu, const double*const rho_s_out, const double*const rho_a_out) const
// {
//     // The (GC) free energy normalized according to the NUMBER OF (ORIGINAL) UNIT CELLS
//     // Each term is already normalized.
//     const double ans = Omega_MF (energies, mu) + mean_V   (rho_s_out, rho_a_out)
//                                                - mean_V_MF(rho_s_out, rho_a_out);
//     return ans;
// }
// 
// double ham1_t::Omega_MF(const double*const energies, const double mu) const
// {
//     // (GC) free energy of the MF Hamiltonian. 
//     // Normalized according to the NUMBER OF (ORIGINAL) UNIT CELLS
//     // It's calculated differently at zero and nonzero temperatures.
//     double accumulator = 0.;
//     
//     if (zerotemp_)
//     { // Remember that there are num_states*num_harmonics energies in kspace.
//       for (int i=0; i<num_states*num_harmonics; ++i)
//         if (energies[i]<=mu)
//           accumulator += (energies[i] - mu);
//     }
//     
//     else // Remember that there are num_states*num_harmonics energies in kspace.
//       for (int i=0; i<num_states*num_harmonics; ++i)
//         accumulator += - T_ * log_1p_exp(-(energies[i] - mu)/T_);
//     
//     // Normalize by the number of unit cells. Also correct for summing over the full BZ.
//     accumulator /= (double)(num_unit_cells*num_harmonics);
//     
//     return accumulator;
// }
// 
// double ham1_t::mean_V(const double*const rho_s_out, const double*const rho_a_out) const
// {
//     // Expectation value of the full interaction term in the Hamiltonian
//     // Normalized according to the NUMBER OF (ORIGINAL) UNIT CELLS
//     
//     // FIRST TERM ***********************************************************************
//     
//     // For convenience, define
//     double rho_A_out [num_harmonics];
//     double rho_B_out [num_harmonics];
//     for (int Q=0; Q<num_harmonics; ++Q) {
//         rho_A_out[Q] = rho_s_out[Q] + rho_a_out[Q];
//         rho_B_out[Q] = rho_s_out[Q] - rho_a_out[Q];
//     }
//     
//     complex<double> term1 = {0.,0.}; // Variable to accumulate on.
//     complex<double> V[states_per_cell*states_per_cell]={0.}; // Declare the array V, to be evaluated for every Q.
//     for (int Q=0; Q<num_harmonics; ++Q) {
//       Assign_V(Q, V); // The matrix V, which depends on Q, gets assigned
//       term1 += 0.5 * (V[0] * std::conj(rho_A_out[Q]) * rho_A_out[Q] + V[1] * std::conj(rho_A_out[Q]) * rho_B_out[Q]
//                     + V[2] * std::conj(rho_B_out[Q]) * rho_A_out[Q] + V[3] * std::conj(rho_B_out[Q]) * rho_B_out[Q]);
//     }
//     
//     if (abs(std::imag(term1))>1.e-15)
//         std::cout << "WARNING: nonzero imaginary part: term1 = " << term1 << std::endl;
//     
//     // SECOND TERM **********************************************************************
//     
//     //complex<double> term2 = {0.,0.}; // Variable to accumulate on.
//     
//     return std::real(term1);// + std::real(term2);
// }
// 
// double ham1_t::mean_V_MF(const double*const rho_s_out, const double*const rho_a_out) const
// {
//     // Expectation value of the M-F interaction term in the Hamiltonian.
//     // Normalized according to the NUMBER OF (ORIGINAL) UNIT CELLS
//     
//     // For convenience, define
//     double rho_A [num_harmonics];
//     double rho_B [num_harmonics];
//     double rho_A_out [num_harmonics];
//     double rho_B_out [num_harmonics];
//     
//     for (int Q=0; Q<num_harmonics; ++Q) // Assign them values from rho_s_ and rho_a_
//     {
//         rho_A[Q]     = rho_s_[Q]    + rho_a_[Q];
//         rho_B[Q]     = rho_s_[Q]    - rho_a_[Q];
//         rho_A_out[Q] = rho_s_out[Q] + rho_a_out[Q];
//         rho_B_out[Q] = rho_s_out[Q] - rho_a_out[Q];
//     }
//     
//     complex<double> ans = {0.,0.}; // Variable for the summation
//     
//     // Declare the array V, to be evaluated for every Q.
//     complex<double> V[states_per_cell*states_per_cell]={0.};
//     for (int Q=0; Q<num_harmonics; ++Q)
//     {
//         Assign_V(Q, V);// The matrix V, which depends on Q, gets assigned
//         ans += V[0] * std::conj(rho_A[Q]) * rho_A_out[Q] + V[1] * std::conj(rho_A[Q]) * rho_B_out[Q]
//              + V[2] * std::conj(rho_B[Q]) * rho_A_out[Q] + V[3] * std::conj(rho_B[Q]) * rho_B_out[Q];
//     }
//     
//     if (abs(std::imag(ans))>1.e-15)
//         std::cout << "WARNING: nonzero imaginary part: ans = " << ans << std::endl;
//     
//     return std::real(ans);
// }
