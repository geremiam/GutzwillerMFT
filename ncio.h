// ncio.h
/* Module containing IO functionalities with NetCDF. */
#ifndef NCIO_H
#define NCIO_H

#include <netcdf.h> // For NetCDF interface

void ErrorHandler(const int status);

class newDS_t
{
private:
    // Private copy constructor (prohibits copy creation)
    newDS_t(const newDS_t&);
    // Private assignment operator (prohibits assignment)
    const newDS_t& operator=(const newDS_t&);
    
    const size_t dims_num_;
    const size_t vars_num_;
    
    int*const coord_varid_; // Array of IDs for coord vars (may not be used in whole or in part)
    int*const varid_; // Array of IDs for user-defined variables
    int varid_energy_=-99; // ID for energy variable
    int varid_mu_=-99; // ID for mu variable (the chemical potential)
    int varid_loops_=-99; // ID for numloops variable
public:
    int ncid_=-99; // ID for dataset
    int*const dimid_; // Array of IDs for dimensions.
    /* The last one (after user-added dims) is complex dim, as has to be the case because 
    complex numbers are stored with real and imaginary parts in consecutive addresses in 
    memory. */
    
    // Constructor declaration
    newDS_t(const size_t dims_num, const std::string*const dim_names, const size_t*const dim_lengths,
            const size_t vars_num, const std::string*const var_names, const bool*const var_complex,
            const std::string GlobalAttr, const std::string path="");
    ~newDS_t(); // Destructor declaration
    
    void DefCoordVar(const int dimindex, const std::string name); // Define a coord variable
    int DefCoordVar_custom(const int dimindex, const std::string name); // Define an individual coord vari whose ID the user keeps track of.
    
    void EndDef(); // Must be called to exit define mode
    
    void WriteCoordVar(const int dimindex, const double*const coord_var); // Write a coordinate variable
    void WriteCoordVar_custom(const int coord_varid, const double*const coord_var); // Write a coord var whose ID the user keeps track of.
    void WriteVars(const double*const*const vars); // Write variables
    void WriteEnergy(const double*const energy); // Write energy variable
    void Writemu(const double*const energy); // Write mu variable
    void WriteLoops(const int*const loops); // Write numloops variable
};

#endif
