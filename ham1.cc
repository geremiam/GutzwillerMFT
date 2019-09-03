// ham1.cc
/* Description */

#include <iostream>
#include <iomanip> // For the function std::setprecision()
#include <sstream> // For stringstreams
#include <complex> // Needed to import alloc.h
#include "alloc.h" // Allocation/deallocation of arrays
#include "array_init.h" // Initialization of arrays
#include "math.h" // Various math functions
#include "misc_routines.h"
#include "chempot.h"
#include "kspace.h" // Defines a class for holding a band structure
#include "diag.h" // Routines for finding evals and evecs
#include "ham1.h" // Include header file for consistency check
//using std::cos; using std::sin; using std::conj; // Not sure if these are necessary
using std::polar;
using std::to_string;
using std::complex;
using std::abs;

double Set_chi(const int counter, const int*const counter_vals, const double*const chi_vals, const int len, const bool with_output=false)
{
    /* When 'counter' reaches one of the values in 'counter_values', the corresponding 
    entry of 'chi_values' is returned. Values in counter_vals must be strictly increasing.*/
    double chi = 1.; // Starting value should be 1.
    for (int i=0; i<len; ++i)
    {
        if (counter>=counter_vals[i]) 
            chi = chi_vals[i];
        if (counter==counter_vals[i]) 
        {
            if (with_output)
                std::cout << "\n\t Counter has reached " << counter
                          << ".\tchi = " << chi << "\n\n";
        }
    }
    
    return chi;
}


int ham1_t::idx(const int Q, const int alpha) const
{
    // Gives the appropriate index in the direct product between Q and alpha.
    // The direct product is defined so that alpha varies faster than Q.
    return idx_composite(num_harmonics, states_per_cell, Q, alpha);
}

void ham1_t::resetMFs()
{
    // Resets MFs to default starting values.
    rho_s_[0] = 0.2;
    rho_a_[0] = 0.1;
    rho_s_[1] = 0.2;
    rho_a_[1] = 0.1;
    rho_s_[2] = 0.0;
    rho_a_[2] = 0.1;
    rho_s_[3] = 0.2;
    rho_a_[3] = 0.1;
}

/* Methods that set parameters with interdependencies. */
void ham1_t::assign_rho(const double rho)
{
    // Assign rho_ and dependent variables. Gets called in the constructor.
    rho_ = rho;
    filled_states = (int)( rho * (double)(num_states) );
}

void ham1_t::set_zerotemp()
{
    zerotemp_ = true; // Make sure this is true.
    T_ = -1.; // Irrelevant when zerotemp_ is true
}

void ham1_t::set_nonzerotemp(const double T)
{
    zerotemp_ = false; // Make sure this is false.
    T_ = T; // Set the private member to the right value
}


/* Constructor and destructor */
ham1_t::ham1_t(const int ka_pts, const int kb_pts, const int kc_pts, const double tol, const int loops_lim,
               const int mixing_vals_len, const int*const counter_vals, const double*const mixing_vals)
    :ka_pts_(ka_pts), kb_pts_(kb_pts), kc_pts_(kc_pts), loops_lim_(loops_lim), tol_(tol)
{
    /* Constructor implementation */
    assign_rho(rho_); // Sets initial value of 'filled_states'
    resetMFs(); // Sets the mean fields to default value
    std::cout << "ham1_t instance created.\n";
    
    const int len_default = 13;
    const int   counter_vals_default [len_default] = { 10,  20,  30,  40,  50,  60,  90, 120, 150,  180,  210,  240,  270};
    const double mixing_vals_default [len_default] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.08, 0.06, 0.04, 0.02};
    
    if (mixing_vals_len==0) // Use the default values
    {
      mixing_vals_len_ = len_default; // The default value
      counter_vals_ = new int [mixing_vals_len_]; // Reserve memory for arrays
      mixing_vals_  = new double [mixing_vals_len_];
      for (int i=0; i<mixing_vals_len_; ++i) // Assign values to default ones
      {
        counter_vals_[i] = counter_vals_default[i];
        mixing_vals_[i]  = mixing_vals_default[i];
      }
    }
    
    else // Use the user input
    {
      mixing_vals_len_ = mixing_vals_len;
      counter_vals_ = new int [mixing_vals_len_]; // Reserve memory for arrays
      mixing_vals_  = new double [mixing_vals_len_];
      for (int i=0; i<mixing_vals_len_; ++i) // Assign values to default ones
      {
        counter_vals_[i] = counter_vals[i];
        mixing_vals_[i]  = mixing_vals[i];
      }
    }
}

ham1_t::~ham1_t()
{
    /* Destructor implementation */
    std::cout << "ham1_t instance deleted.\n";
}

// Methods for assigning the Hamiltonian
void ham1_t::Assign_h(double ka, double kb, double kc, complex<double>*const h) const
{
    // With this convention, it is not necessary to use momenta inside the FBZ.
    //ka = FBZ(ka, a_);
    //kb = FBZ(kb, a_);
    //kc = FBZ(kc, c_);
    // Assign to h
    h[0] = -2.*t1p_ * (cos(a_*ka) + cos(a_*kb)) - 2.*t3_ * cos(c_*kc);
    h[1] = //(- t1_ * zeta(ka,kb) - t2A_ * f(ka,kb,kc) - t2B_ * std::conj(f(ka,kb,kc))) * polar(1.,-a_*(ka+kb)/2.);
           - t1_  * ( 1. + polar(1.,-a_*ka) + polar(1.,-a_*kb) + polar(1.,-a_*ka - a_*kb) )
           - t2A_ * ( polar(1.,-c_*kc) + polar(1.,+c_*kc - a_*ka) + polar(1.,+c_*kc - a_*kb) + polar(1.,-a_*ka - a_*kb - c_*kc) )
           - t2B_ * ( polar(1.,+c_*kc) + polar(1.,-c_*kc - a_*ka) + polar(1.,-c_*kc - a_*kb) + polar(1.,-a_*ka - a_*kb + c_*kc) );
    h[2] = std::conj(h[1]);
    h[3] = h[0];
}
void ham1_t::Assign_V(const int Q, complex<double>*const V) const
{
    // Assigns the matrix to V with qa = Qa[Q], qb = Qb[Q], etc.
    const double qa = Qa[Q];
    const double qb = Qb[Q];
    const double qc = Qc[Q];
    // Assign to V
    V[0] = 2.*V1p_ * (cos(a_*qa) + cos(a_*qb)) + 2.*V3_ * cos(c_*qc);
    V[1] = (V1_ + 2.*V2_*cos(c_*qc)) * (1. + polar(1.,-a_*qa) + polar(1.,-a_*qb) + polar(1.,-a_*qa-a_*qb));
    V[2] = std::conj(V[1]);
    V[3] = V[0];
}
/*void ham1_t::Assign_V_manual(const int Q, complex<double>*const V) const
{
    // Manually computed assignment formulas for V. Assigns the matrix to V with 
    // qa = Qa[Q], qb = Qb[Q], etc.
    // NEEDS TO BE REWRITTEN FOR DIFFERENT SETS OF HARMONICS
    if (Q==0) {
      V[0] = 4.*V1p_ + 2.*V3_;
      V[1] = 4.*(V1_ + 2.*V2_);
    }
    
    else if (Q==1) {
      V[0] = 2.*V3_;
      V[1] = 0.;
    }
    
    else if (Q==2) {
      V[0] = 2.*V3_;
      V[1] = 0.;
    }
    
    else if (Q==3) {
      V[0] = -4.*V1p_ + 2.*V3_;
      V[1] = 0.;
    }
    else
      std::cout << "WARNING: In routine Assign_V_manual, index Q is out of range." << std::endl;
    
    V[2] = std::conj(V[1]);
    V[3] = V[0];
}*/

// Assigns Hamiltonian
void ham1_t::Assign_ham(const double ka, const double kb, const double kc, complex<double>*const*const ham_array) const
{
    /* Given the arguments momentum as well as member parameters, calculate the 2*2
    k-space Hamiltonian and assign it to ham_array in full storage layout. The diag 
    routine only uses the lower triangle, but we assign the whole thing. */
    
    complex<double>*const*const& H = ham_array; // For convenience, define reference
    // Initialize elements to zero
    for (int i=0; i<num_bands*num_bands; ++i)
        H[0][i] = 0.;
    // **********************************************************************************
    
    // Add the kinetic terms
    
    // Array for holding the 2*2 kinetic part of the Hamiltonian
    complex<double> h[states_per_cell*states_per_cell]={0.};
    for (int Q=0; Q<num_harmonics; ++Q)
    {
      Assign_h(ka+Qa[Q], kb+Qb[Q], kc+Qc[Q], h); // Assign the 2*2 kinetic part of the Hamiltonian
      H[idx(Q,0)][idx(Q,0)] += h[0];
      H[idx(Q,0)][idx(Q,1)] += h[1];
      H[idx(Q,1)][idx(Q,0)] += h[2];
      H[idx(Q,1)][idx(Q,1)] += h[3];
    }
    // **********************************************************************************
    
    
    // Declare an array rhotilde for each alpha. These have one index for the Harmonic.
    complex<double> rhotilde_A [num_harmonics]={0.};
    complex<double> rhotilde_B [num_harmonics]={0.};
    // Declare the array V, to be evaluated for every Q, needed to compte rhotilde.
    complex<double> V[states_per_cell*states_per_cell]={0.};
    for (int Q=0; Q<num_harmonics; ++Q)
    {
        // For convenience, assign 
        const double rho_A = rho_s_[Q] + rho_a_ [Q];
        const double rho_B = rho_s_[Q] - rho_a_ [Q];
        Assign_V(Q, V);// The matrix V, which depends on Q, gets assigned
        rhotilde_A[Q] = V[0]*rho_A + V[1]*rho_B;// V[0,0]*rho_A_[Q] + V[0,1]*rho_B_[Q];
        rhotilde_B[Q] = V[2]*rho_A + V[3]*rho_B;// V[1,0]*rho_A_[Q] + V[1,1]*rho_B_[Q];
    }
    
    // Add the density terms to the Hamiltonian
    for (int P=0; P<num_harmonics; ++P)
      for (int Q=0; Q<num_harmonics; ++Q)
      {
          const int Q_sum = addition_table[P][Q];
          // Loop over the momentum double sum; use the addition table.
          H[idx(P,0)][idx(Q_sum,0)] += std::conj(rhotilde_A[Q]);
          H[idx(P,1)][idx(Q_sum,1)] += std::conj(rhotilde_B[Q]);
      }
}


// Methods for computing mean fields
void ham1_t::AddContribution_rho(const double*const occs, const complex<double>*const*const evecs, double*const rho_A, double*const rho_B) const
{
    // Evaluates the contribution to rho_A and rho_B from a single k (see notes) and adds 
    // it to the arrays "rho_A" and "rho_B".
    
    for (int Q=0; Q<num_harmonics; ++Q) // Loop over the harmonics i.e. the components of rho_A and rho_B
    {
      complex<double> temp_A = {0.,0.}; // These will be used in the summation.
      complex<double> temp_B = {0.,0.}; // Initialize to zero
      
      for (int beta=0; beta<states_per_cell; ++beta) // summation loop
        for (int P=0; P<num_harmonics; ++P)          // summation loop
          for (int R=0; R<num_harmonics; ++R)        // summation loop
          {
            const int QplusR = addition_table[Q][R];
            
            temp_A += conj(evecs[idx(R,0)][idx(P,beta)]) * evecs[idx(QplusR,0)][idx(P,beta)] * occs[idx(P,beta)];
            temp_B += conj(evecs[idx(R,1)][idx(P,beta)]) * evecs[idx(QplusR,1)][idx(P,beta)] * occs[idx(P,beta)];
          }
      if (abs(std::imag(temp_A))>1.e-15)
        std::cout << "WARNING: nonzero imaginary part: temp_A = " << temp_A << std::endl;
      if (abs(std::imag(temp_B))>1.e-15)
        std::cout << "WARNING: nonzero imaginary part: temp_B = " << temp_B << std::endl;
      // Normalize by the number of unit cells, plus an extra factor of num_harmonics 
      // accounting for the fact that we are summing over the FULL BZ instead of the RBZ.
      rho_A[Q] += std::real(temp_A) / (double)(num_unit_cells*num_harmonics);
      rho_B[Q] += std::real(temp_B) / (double)(num_unit_cells*num_harmonics);
    }
}


void ham1_t::ComputeMFs_old(double*const rho_s_out, double*const rho_a_out, double*const HFE_p, double*const mu_p) const
{
    // Declare (and construct) an instance of kspace_t.
    // *** The sum can be over the full BZ (instead of the RBZ) as long as this is 
    // compensated for in the calculation of the MFS and of the chemical potential. ***
    kspace_t kspace(a_, a_, c_, ka_pts_, kb_pts_, kc_pts_, num_bands);
    
    // Step 1: diagonalize to find all the energy evals and store them in kspace
    #pragma omp parallel default(none) shared(kspace)
    {
    // array to hold Ham (local to thread)
    complex<double>*const*const ham_array = Alloc2D_z(ham_array_rows, ham_array_cols);
    ValInitArray(ham_array_rows*ham_array_cols, &(ham_array[0][0])); //Initialize to zero
    #pragma omp for collapse(3)
    // Given the parameters, diagonalize the Hamiltonian at each grid point
    for (int i=0; i<ka_pts_; ++i)
      for (int j=0; j<kb_pts_; ++j)
        for (int k=0; k<kc_pts_; ++k)
        {
          Assign_ham(kspace.ka_grid[i], kspace.kb_grid[j], kspace.kc_grid[k], ham_array);
          const int index = kspace.index(i, j, k, 0); // Last index is the band index
          simple_zheev(num_bands, &(ham_array[0][0]), &(kspace.energies[index]));
        }
    Dealloc2D(ham_array);
    }
    
    // Step 2: Use all energies to compute chemical potential
    // WE MUST CORRECT FOR THE OVERCOUNTED BZ INTEGRATION BY INCREASING num_states AND 
    // filled_states BY A FACTOR OF num_harmonics. OTHERWISE, ONLY PART OF THE kspace IS
    // USED TO CALCULATE mu.
    double mu = 666.;
    if (zerotemp_) // (ELEMENTS GET REORDERED)
        mu = FermiEnerg(num_states*num_harmonics, filled_states*num_harmonics, kspace.energies);
    else // USES OPENMP PARALLELIZATION
    {
        const bool show_output = false;
        const bool usethreads = true;
        mu = ChemPotBisec(num_states*num_harmonics, filled_states*num_harmonics, kspace.energies, T_, show_output, usethreads);
    }
    
    
    // Step 3: Use all the occupation numbers and the evecs to find the order parameter
    // Probably best to diagonalize a second time to avoid storing the evecs
    double rho_A_accum [num_harmonics] = {0.}; // IMPORTANT: MUST BE INITIALIZED TO ZERO
    double rho_B_accum [num_harmonics] = {0.}; // IMPORTANT: MUST BE INITIALIZED TO ZERO
    
    #pragma omp parallel default(none) firstprivate(mu) shared(kspace) reduction(+:rho_A_accum[0:num_harmonics], rho_B_accum[0:num_harmonics])
    {
    complex<double>*const*const ham_array = Alloc2D_z(ham_array_rows, ham_array_cols); // array to hold Ham (local to thread)
    complex<double>*const*const     evecs = Alloc2D_z(num_bands, num_bands); // Array to hold evecs (local to each thread)
    double*const evals = new double [num_bands]; // Array to hold evals in second loop (local to thread)
    double*const  occs = new double [num_bands]; // Array to hold occupations (local to thread)
    ValInitArray(ham_array_rows*ham_array_cols, &(ham_array[0][0])); // Initialize to zero
    ValInitArray(num_bands*num_bands, &(evecs[0][0])); // Initialize to zero
    ValInitArray(num_bands, evals); // Initialize to zero
    ValInitArray(num_bands, occs); // Initialize to zero
    
    #pragma omp for collapse(3)
    for (int i=0; i<ka_pts_; ++i)
      for (int j=0; j<kb_pts_; ++j)
        for (int k=0; k<kc_pts_; ++k)
        {
          Assign_ham(kspace.ka_grid[i], kspace.kb_grid[j], kspace.kc_grid[k], ham_array);
          simple_zheev(num_bands, &(ham_array[0][0]), evals, true, &(evecs[0][0]));
          // Calculate occupations from energies, mu, and temperature
          Occupations(num_bands, mu, evals, occs, zerotemp_, T_);
          
          AddContribution_rho(occs, evecs, rho_A_accum, rho_B_accum);
        }
    delete [] occs; // Deallocate memory for arrays. 
    delete [] evals;
    Dealloc2D(evecs);
    Dealloc2D(ham_array);
    }
    
    // Assign the results to the output arrays
    
    for (int Q=0; Q<num_harmonics; ++Q)
    {
        rho_s_out[Q] = 0.5 * (rho_A_accum[Q] + rho_B_accum[Q]);
        rho_a_out[Q] = 0.5 * (rho_A_accum[Q] - rho_B_accum[Q]);
    }
    
    // The kspace_t destructor is called automatically.
    
    if (HFE_p!=NULL)
      *HFE_p = Helmholtz(kspace.energies, mu, rho_s_out, rho_a_out);
    if (mu_p!=NULL)
      *mu_p  = mu;
}

void ham1_t::ComputeMFs    (double*const rho_s_out, double*const rho_a_out, double*const HFE_p, double*const mu_p) const
{
    // Declare (and construct) an instance of kspace_t.
    // *** The sum can be over the full BZ (instead of the RBZ) as long as this is 
    // compensated for in the calculation of the MFS and of the chemical potential. ***
    const bool with_output=false; const bool with_evecs=true;
    kspace_t kspace(a_, a_, c_, ka_pts_, kb_pts_, kc_pts_, num_bands, with_output, with_evecs);
    
    // Step 1: diagonalize to find all the energy evals and evecs and store them in kspace
    #pragma omp parallel default(none) shared(kspace)
    {
    // array to hold Ham (local to thread)
    complex<double>*const*const ham_array = Alloc2D_z(ham_array_rows, ham_array_cols);
    ValInitArray(ham_array_rows*ham_array_cols, &(ham_array[0][0])); //Initialize to zero
    #pragma omp for collapse(3)
    // Given the parameters, diagonalize the Hamiltonian at each grid point
    for (int i=0; i<ka_pts_; ++i)
      for (int j=0; j<kb_pts_; ++j)
        for (int k=0; k<kc_pts_; ++k)
        {
          Assign_ham(kspace.ka_grid[i], kspace.kb_grid[j], kspace.kc_grid[k], ham_array);
          const int evals_ind = kspace.index(i, j, k, 0); // Last index is the band index
          const int     k_ind = kspace.k_ind(i, j, k);
          simple_zheev(num_bands, &(ham_array[0][0]), &(kspace.energies[evals_ind]), true, &(kspace.evecs[k_ind][0][0]));
        }
    Dealloc2D(ham_array);
    }
    
    // Step 2: Use all energies to compute chemical potential
    // WE MUST CORRECT FOR THE OVERCOUNTED BZ INTEGRATION BY INCREASING num_states AND 
    // filled_states BY A FACTOR OF num_harmonics. OTHERWISE, ONLY PART OF THE kspace IS
    // USED TO CALCULATE mu.
    double mu = 666.;
    if (zerotemp_) // Array kspace.energies is left unchanged.
        mu = FermiEnerg_cpy(num_states*num_harmonics, filled_states*num_harmonics, kspace.energies);
    else // USES OPENMP PARALLELIZATION
    {
        const bool show_output = false;
        const bool usethreads = true;
        mu = ChemPotBisec(num_states*num_harmonics, filled_states*num_harmonics, kspace.energies, T_, show_output, usethreads);
    }
    
    
    // Step 3: Use all the occupation numbers and the evecs to find the order parameter
    // Probably best to diagonalize a second time to avoid storing the evecs
    double rho_A_accum [num_harmonics] = {0.}; // IMPORTANT: MUST BE INITIALIZED TO ZERO
    double rho_B_accum [num_harmonics] = {0.}; // IMPORTANT: MUST BE INITIALIZED TO ZERO
    
    // https://www.nersc.gov/assets/OMP-common-core-SC17.pdf for reduction on array sections
    // https://www.archer.ac.uk/training/virtual/2017-09-27-OpenMP4/OpenMP45VT.pdf
    #pragma omp parallel default(none) firstprivate(mu) shared(kspace) reduction(+:rho_A_accum[0:num_harmonics], rho_B_accum[0:num_harmonics])
    {
    double*const  occs = new double [num_bands]; // Array to hold occupations (local to thread)
    ValInitArray(num_bands, occs); // Initialize to zero
    
    #pragma omp for collapse(3)
    for (int i=0; i<ka_pts_; ++i)
      for (int j=0; j<kb_pts_; ++j)
        for (int k=0; k<kc_pts_; ++k)
        {
          // Calculate occupations from energies, mu, and temperature
          const int evals_ind = kspace.index(i, j, k, 0); // Last index is the band index
          const int     k_ind = kspace.k_ind(i, j, k);
          Occupations(num_bands, mu, &(kspace.energies[evals_ind]), occs, zerotemp_, T_);
          AddContribution_rho(occs, kspace.evecs[k_ind], rho_A_accum, rho_B_accum);
        }
    delete [] occs; // Deallocate memory for arrays. 
    }
    
    // Assign the results to the output arrays
    for (int Q=0; Q<num_harmonics; ++Q)
    {
        rho_s_out[Q] = 0.5 * (rho_A_accum[Q] + rho_B_accum[Q]);
        rho_a_out[Q] = 0.5 * (rho_A_accum[Q] - rho_B_accum[Q]);
    }
    
    // The kspace_t destructor is called automatically.
    if (HFE_p!=NULL) // The free energy is assigned to HFE_p
      *HFE_p = Helmholtz(kspace.energies, mu, rho_s_out, rho_a_out);
    if (mu_p!=NULL) // The chemical potential is assigned to mu_p
      *mu_p  = mu;
}



std::string ham1_t::GetAttributes()
{
    // Define a string of metadata
    // Watch out: this method resets the values of the MFs to their starting values.
    resetMFs();
        
    std::ostringstream strs; // Declare a stringstream to which to write the attributes
    strs << "Gyrotropic SSB (ham1)"
         << "\na = " << a_ << ", c = " << c_
         << "\nk-space points: " << ka_pts_ << "," << kb_pts_ << "," << kc_pts_ << "; num_bands = " << num_bands
         << "\nzerotemp = " << zerotemp_ << ", T = " << T_ 
         << "\ntol = " << tol_ << ", loops_lim = " << loops_lim_
         << "\nrho = " << rho_
         << "\nt1 = " << t1_ << ", t1p = " << t1p_ << ", t2A = " << t2A_ << ", t2B = " << t2B_ << ", t3 = " << t3_ 
         << "\nV1 = " << V1_ << "; V1p = " << V1p_ << ", V2 = " << V2_ << ", V3 = " << V3_;
    
    strs << "\nList of harmonics: ";
    for (int Q=0; Q<num_harmonics; ++Q)
        strs << std::endl << "Q" << Q << " = (" << Qa[Q] << ", " << Qb[Q] << ", " << Qc[Q] << ")";
    
    strs << "\nStarting values of the MFs: ";
    for (int Q=0; Q<num_harmonics; ++Q)
        strs << std::endl << "rho_s[" << Q << "] = " << rho_s_[Q] << ", rho_a[" << Q << "] = " << rho_a_[Q];
        
    const std::string Attributes = strs.str(); // Write the stream contents to a string
    
    return Attributes;
}


bool ham1_t::FixedPoint(int*const num_loops_p, const bool with_output)
{
    // Performs the iterative self-consistent search using the parameters from ham1. 
    // The initial values of the MF attributes are used as the starting values for the 
    // search; the end values are also stored in these same MF attributes.
    if (with_output) // Format display output and set precision
        std::cout << std::scientific << std::showpos;
    
    // Manually set rho_s_[0] to rho_, so the starting value for this "MF" is correct.
    rho_s_[0] = rho_;
    
    // Declare output variables and INITIALIZE THEM TO INPUT VALUES
    double rho_s_out [num_harmonics];
    double rho_a_out [num_harmonics];
    copy_array(num_harmonics, rho_s_, rho_s_out); // Copy rho_s_ into rho_s_out
    copy_array(num_harmonics, rho_a_, rho_a_out); // Copy rho_a_ into rho_a_out
    double HFE_prev = 0.; // For keeping track of free energy at previous step
    
    if (with_output)
      std::cout << "\t" "rho_s_" "\t\t" "rho_a_" << std::endl;
    
    int counter = 0; // Define counter for number of loops
    bool converged=false, fail=false; // Used to stop the while looping
    do // Iterate until self-consistency is achieved
    {
        ++counter; // Increment counter
        
        
        // Past a certain number of loops, we mix in part of the previous input vals
        // Mixing fraction chi (chi=1 corresponds to using fully new value)
        const double chi = Set_chi(counter, counter_vals_, mixing_vals_, mixing_vals_len_, with_output);
        
        for (int Q=0; Q<num_harmonics; ++Q) // Update mean-field values
        {
            rho_s_[Q] = (1.-chi)*rho_s_[Q] + chi*rho_s_out[Q];
            rho_a_[Q] = (1.-chi)*rho_a_[Q] + chi*rho_a_out[Q];
        }
        HFE_prev = HFE_; // Store previous free energy in HFE_prev
        
        
        if (with_output)
        {
          std::cout << counter << "  chi = " << chi << "\ninput\t";
          for (int Q=0; Q<num_harmonics; ++Q)
          {
            std::cout << rho_s_[Q] << "\t";
            std::cout << rho_a_[Q] << "\t";
          }
        }
        
        // Compute MFs. Output is assigned to the 'out' arguments.
        // The HFE for the final set of MF values is assigned to the attribute HFE_.
        // Likewise, the final chemical potential is assigned to mu_.
        ComputeMFs(rho_s_out, rho_a_out, &HFE_, &mu_);
        
        if (with_output)
          std::cout << "|  " << HFE_ << std::endl;
        
        double rho_s_diff [num_harmonics]; // Declare arrays to store the changes in MF
        double rho_a_diff [num_harmonics];
        for (int Q=0; Q<num_harmonics; ++Q) // Differences between outputs and inputs
        {
            rho_s_diff[Q] = rho_s_out[Q] - rho_s_[Q];
            rho_a_diff[Q] = rho_a_out[Q] - rho_a_[Q];
        }
        
        if (with_output) // Print the differences
        {
          std::cout << "diff\t";
          for (int Q=0; Q<num_harmonics; ++Q)
          {
            std::cout << rho_s_diff[Q] << "\t";
            std::cout << rho_a_diff[Q] << "\t";
          }
          std::cout << "|  " << HFE_-HFE_prev << std::endl;
        }
        
        // Test for convergence. The density (first element of rho_s_diff) is ignored.
        converged = check_bound_array(num_harmonics-1, tol_, rho_s_diff+1) 
                 && check_bound_array(num_harmonics, tol_, rho_a_diff);
        fail = (!converged) && (counter>loops_lim_); // Must come after converged line
        
    } while (!converged && !fail);
    
    // Unless num_loops_p is the null pointer, assign the number of loops to its location
    if (num_loops_p!=NULL) *num_loops_p = counter;
    
    // We make sure that either converged or fail is true.
    if ((converged==true) && (fail==true) )
        std::cout << "OUPS 1: This option shouldn't have occurred! (A)\n";
    if ((converged==false) && (fail==false) )
        std::cout << "OUPS 1: This option shouldn't have occurred! (B)\n";
    
    return fail;
}


// **************************************************************************************
// ROUTINES FOR CALCULATING THE FREE ENERGY

double ham1_t::Helmholtz(const double*const energies, const double mu, const double*const rho_s_out, const double*const rho_a_out) const
{
    // The Helmholtz FE normalized according to the NUMBER OF (ORIGINAL) UNIT CELLS
    // Should be used to compare states of a fixed particle number
    // Second term is + mu n (properly normalized).
    // We use the "real" density rho_s_[0] instead of the target density rho_. These may 
    // be different because of filled_states gets rounded to an int.
    const double ans = Omega_trial(energies,mu,rho_s_out,rho_a_out) + mu*rho_s_[0]*states_per_cell;
    return ans;
}

double ham1_t::Omega_trial(const double*const energies, const double mu, const double*const rho_s_out, const double*const rho_a_out) const
{
    // The (GC) free energy normalized according to the NUMBER OF (ORIGINAL) UNIT CELLS
    // Each term is already normalized.
    const double ans = Omega_MF (energies, mu) + mean_V   (rho_s_out, rho_a_out)
                                               - mean_V_MF(rho_s_out, rho_a_out);
    return ans;
}

double ham1_t::Omega_MF(const double*const energies, const double mu) const
{
    // (GC) free energy of the MF Hamiltonian. 
    // Normalized according to the NUMBER OF (ORIGINAL) UNIT CELLS
    // It's calculated differently at zero and nonzero temperatures.
    double accumulator = 0.;
    
    if (zerotemp_)
    { // Remember that there are num_states*num_harmonics energies in kspace.
      for (int i=0; i<num_states*num_harmonics; ++i)
        if (energies[i]<=mu)
          accumulator += (energies[i] - mu);
    }
    
    else // Remember that there are num_states*num_harmonics energies in kspace.
      for (int i=0; i<num_states*num_harmonics; ++i)
        accumulator += - T_ * log_1p_exp(-(energies[i] - mu)/T_);
    
    // Normalize by the number of unit cells. Also correct for summing over the full BZ.
    accumulator /= (double)(num_unit_cells*num_harmonics);
    
    return accumulator;
}

double ham1_t::mean_V(const double*const rho_s_out, const double*const rho_a_out) const
{
    // Expectation value of the full interaction term in the Hamiltonian
    // Normalized according to the NUMBER OF (ORIGINAL) UNIT CELLS
    
    // FIRST TERM ***********************************************************************
    
    // For convenience, define
    double rho_A_out [num_harmonics];
    double rho_B_out [num_harmonics];
    for (int Q=0; Q<num_harmonics; ++Q) {
        rho_A_out[Q] = rho_s_out[Q] + rho_a_out[Q];
        rho_B_out[Q] = rho_s_out[Q] - rho_a_out[Q];
    }
    
    complex<double> term1 = {0.,0.}; // Variable to accumulate on.
    complex<double> V[states_per_cell*states_per_cell]={0.}; // Declare the array V, to be evaluated for every Q.
    for (int Q=0; Q<num_harmonics; ++Q) {
      Assign_V(Q, V); // The matrix V, which depends on Q, gets assigned
      term1 += 0.5 * (V[0] * std::conj(rho_A_out[Q]) * rho_A_out[Q] + V[1] * std::conj(rho_A_out[Q]) * rho_B_out[Q]
                    + V[2] * std::conj(rho_B_out[Q]) * rho_A_out[Q] + V[3] * std::conj(rho_B_out[Q]) * rho_B_out[Q]);
    }
    
    if (abs(std::imag(term1))>1.e-15)
        std::cout << "WARNING: nonzero imaginary part: term1 = " << term1 << std::endl;
    
    // SECOND TERM **********************************************************************
    
    //complex<double> term2 = {0.,0.}; // Variable to accumulate on.
    
    return std::real(term1);// + std::real(term2);
}

double ham1_t::mean_V_MF(const double*const rho_s_out, const double*const rho_a_out) const
{
    // Expectation value of the M-F interaction term in the Hamiltonian.
    // Normalized according to the NUMBER OF (ORIGINAL) UNIT CELLS
    
    // For convenience, define
    double rho_A [num_harmonics];
    double rho_B [num_harmonics];
    double rho_A_out [num_harmonics];
    double rho_B_out [num_harmonics];
    
    for (int Q=0; Q<num_harmonics; ++Q) // Assign them values from rho_s_ and rho_a_
    {
        rho_A[Q]     = rho_s_[Q]    + rho_a_[Q];
        rho_B[Q]     = rho_s_[Q]    - rho_a_[Q];
        rho_A_out[Q] = rho_s_out[Q] + rho_a_out[Q];
        rho_B_out[Q] = rho_s_out[Q] - rho_a_out[Q];
    }
    
    complex<double> ans = {0.,0.}; // Variable for the summation
    
    // Declare the array V, to be evaluated for every Q.
    complex<double> V[states_per_cell*states_per_cell]={0.};
    for (int Q=0; Q<num_harmonics; ++Q)
    {
        Assign_V(Q, V);// The matrix V, which depends on Q, gets assigned
        ans += V[0] * std::conj(rho_A[Q]) * rho_A_out[Q] + V[1] * std::conj(rho_A[Q]) * rho_B_out[Q]
             + V[2] * std::conj(rho_B[Q]) * rho_A_out[Q] + V[3] * std::conj(rho_B[Q]) * rho_B_out[Q];
    }
    
    if (abs(std::imag(ans))>1.e-15)
        std::cout << "WARNING: nonzero imaginary part: ans = " << ans << std::endl;
    
    return std::real(ans);
}
