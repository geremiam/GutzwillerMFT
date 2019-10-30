// ncio.cc

#include <iostream> // For IO to command line
#include <string> // c++ strings
#include <time.h> // time_t, time, ctime
#include "ncio.h" // Include header file for consistency check

// Documentation: https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-c/

// Internal routine declarations
std::string DateTime();
// ***************
// ***************

newDS_t::newDS_t(const size_t dims_num, const std::string*const dim_names, const size_t*const dim_lengths,
                 const size_t vars_num, const std::string*const var_names, const bool*const var_complex,
                 const std::string GlobalAttr, const std::string path, const bool prune)
    :dims_num_(dims_num), dimid_(new int [dims_num+1]), coord_varid_(new int [dims_num]), 
     vars_num_(vars_num), varid_(new int [vars_num]), 
     prune_(prune)
{
    std::cout << "newDS_t instance created.\n";
    
    // Creation of the dataset
    // Existing files of the same name are overwritten. We avoid using NetCDF4 for compatibility with the python interface.
    const std::string Filename = path + DateTime();
    ErrorHandler( nc_create(Filename.c_str(), NC_CLOBBER, &ncid_) ); // The dataset ID is assigned to ncid. 
    
    // Global attribute
    // The string GlobalAttr is set as a global attribute of the dataset.
    const size_t len = GlobalAttr.length() + 1; // Length of char array (with terminator)
    ErrorHandler( nc_put_att_text(ncid_, NC_GLOBAL, "GlobalAttributes", len, GlobalAttr.c_str()) );
    
    // Creation of dimensions
    // Their names and lengths are supplied as arguments.
    for (int i=0; i<dims_num_; ++i)
        ErrorHandler( nc_def_dim(ncid_, dim_names[i].c_str(), dim_lengths[i], &(dimid_[i])) );
    
    // Extra dimension for real and imaginary parts (of length 2).
    // Must be last component of dimid_, because this is how complex numbers are stored.
    ErrorHandler( nc_def_dim(ncid_, "complex", 2, &(dimid_[dims_num_])) );
    
    // Declaration of data variables
    // We presume they each span all dimensions. Complex vars also span the dimension 
    // "complex" as their last dimension. This has to be the last dimension because complex 
    // numbers are stored with real and imaginary parts in consecutive memory addresses.
    for (int j=0; j<vars_num_; ++j)
    {
        if (var_complex[j]==false) // In this case the variable is real
          ErrorHandler( nc_def_var(ncid_, var_names[j].c_str(), NC_DOUBLE, dims_num_,   dimid_, &(varid_[j])) );
        else // In this case the variable is complex
          ErrorHandler( nc_def_var(ncid_, var_names[j].c_str(), NC_DOUBLE, dims_num_+1, dimid_, &(varid_[j])) );
    }
    
    // Variable for the number of loops
    ErrorHandler( nc_def_var(ncid_, "numloops", NC_INT, dims_num_, dimid_, &varid_loops_ ));
    
}

void newDS_t::DefCoordVar(const int dimindex, const std::string name)
{
    // Definition of individual coord variables
    ErrorHandler( nc_def_var(ncid_, name.c_str(), NC_DOUBLE, 1, &(dimid_[dimindex]), &(coord_varid_[dimindex])) );
}

int newDS_t::DefCoordVar_custom(const int dimindex, const std::string name)
{
    // Definition of individual coord variables whose ID the user keeps track of.
    int coord_varid=-99;
    ErrorHandler( nc_def_var(ncid_, name.c_str(), NC_DOUBLE, 1, &(dimid_[dimindex]), &coord_varid) );
    return coord_varid;
}

void newDS_t::EndDef()
{
    // Exit define mode and enter data mode.
    ErrorHandler( nc_enddef(ncid_) );
}


void newDS_t::WriteCoordVar(const int dimindex, const double*const coord_var)
{
    // Write a single coordinate. The NetCDF C interface expects row-major layout.
    ErrorHandler( nc_put_var_double(ncid_, coord_varid_[dimindex], coord_var) );
}
void newDS_t::WriteCoordVar_custom(const int coord_varid, const double*const coord_var)
{
    // Write a single coordinate. The NetCDF C interface expects row-major layout.
    ErrorHandler( nc_put_var_double(ncid_, coord_varid, coord_var) );
}
void newDS_t::WriteVars(const double*const*const vars)
{
    // Write variables. The NetCDF C interface expects row-major layout.
    // Convert pointers to complex arrays to double* (pointing to first real part).
    for (int j=0; j<vars_num_; ++j)
        ErrorHandler( nc_put_var_double(ncid_, varid_[j], vars[j]) );
}
void newDS_t::WriteLoops(const int*const loops)
{
    // Write energy variable. The NetCDF C interface expects row-major layout.
    ErrorHandler( nc_put_var_int(ncid_, varid_loops_, loops) );
}

newDS_t::~newDS_t()
{
    std::cout << "newDS_t instance deleted.\n";
    
    /* Closes the datased referred to by ncid. */
    ErrorHandler( nc_close(ncid_) );
    
    // Deallocate memory
    delete [] dimid_;
    delete [] coord_varid_;
    delete [] varid_;
}


// ***************

// Internal routines
void ErrorHandler(const int status)
{
    /*
        INTERNAL ROUTINE
        Handles NetCDF errors using the function nc_strerror() from the NetCDF C 
        interface. If the value of retval indicates an error, a message is output.
    */
    if (status != NC_NOERR)
    {
        std::cerr << "\n***NetCDF error***\n"
                  << nc_strerror(status) << "\n\n";
    }
}

std::string DateTime()
{
    /* Returns a c++ string with the current local time and date in the specified format.
    http://www.cplusplus.com/reference/ctime/strftime/ */
    time_t rawtime;
    struct tm* timeinfo;
    char buffer [100];
    
    time(&rawtime); // Assign time to rawtime
    timeinfo = localtime(&rawtime); // Use this to assign timeinfo
    
    // Use this to assign the time with the proper format to buffer
    strftime(buffer,100,"%Y-%m-%d_%Hh%Mmin%Ss",timeinfo);
    
    // Define a c++ string from buffer
    std::string TimeString = buffer;
    
    return TimeString;
}
