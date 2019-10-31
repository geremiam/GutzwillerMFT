// ncio.h

// Module for creating a netCDF dataset in which all variables are defined on the same 
// dimensions. Complex variables are dealt with by using an additional dimension called 
// "complex". Pruning functionality allows smart handling of dimensions of length 1.

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
    
    int ncid_=-99; // ID for dataset
    int*const dimid_; // Array of IDs for dimensions.
    int*const dimid_relevant_; // Array of "relevant" dimension IDs. Different from dimid_ if prune_==true and there are length-1 dimensions.
    // The last one (after user-added dims) is complex dim, as has to be the case because 
    // complex numbers are stored with real and imaginary parts in consecutive addresses in 
    // memory.
    
    const size_t dims_num_; // Number of dimensions
          size_t dims_num_relevant_=-1; // Number of relevant dimensions. Different from dims_num_ if prune==true and there are length-1 dimensions.
          size_t*const dim_lengths_; // To be copied from input
    const size_t vars_num_; // Number of variables
    
    
    const bool prune_; // Activates pruning of length-1 dimensions
    
    // (A given element of coord_varid_ is used only if the user supplies a coord var for that dimension)
    int*const coord_varid_; // Array of IDs for coord vars. Length dims_num_.
    int*const varid_; // Array of IDs for user-defined variables. Length vars_num_.
    int varid_loops_=-99; // ID for numloops variable
    
    void find_relevant_dimensions(); // Used to assign values to dimid_relevant_ from dimid_.
    
public:
    
    // Constructor declaration
    newDS_t(const size_t dims_num, const std::string*const dim_names, const size_t*const dim_lengths,
            const size_t vars_num, const std::string*const var_names, const bool*const var_complex,
            const std::string GlobalAttr, const std::string path="", const bool prune=false);
    ~newDS_t(); // Destructor declaration
    
    void DefCoordVar(const int dimindex, const std::string name); // Define a coord variable
    int  DefCoordVar_custom(const int dimindex, const std::string name); // Define an individual coord var whose ID the user keeps track of.
    
    void EndDef(); // Must be called to exit define mode
    
    void WriteCoordVar(const int dimindex, const double*const coord_var); // Write a coordinate variable
    void WriteCoordVar_custom(const int coord_varid, const double*const coord_var); // Write a coord var whose ID the user keeps track of.
    void WriteVars(const double*const*const vars); // Write variables
    void WriteLoops(const int*const loops); // Write numloops variable
};

#endif
