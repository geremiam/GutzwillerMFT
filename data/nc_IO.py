# nc_IO.py
""" Routines for interacting with NetCDF datasets. """

import numpy as np
from scipy.io import netcdf
from collections import OrderedDict
from sys import argv

def nc_read_var(filename, varname):
    """ Reads a simple format of NetCDF datasets with n dimensions (each having a 
    coordinate variable) and m (non-coordinate) variables. The dimension names and 
    lengths, the coordinate variables, the variable names and the variables are returned. 
    Note that variables are loaded into memory and returned as arrays, so this is not 
    adequate for very large datasets. """
    np.set_printoptions(linewidth=1000)
    
    with netcdf.netcdf_file(filename, 'r') as f: # Open the NetCDF file
        print("Dimensions of the dataset: {}".format(f.dimensions)) # Print all dims
        
        dataset_var = f.variables[varname] # Get the desired variable
        
        var_dims = dataset_var.dimensions # Tuple of dimensions this var is defined on
        print("Dimensions for this variable: {}".format(var_dims))
        
        var = dataset_var.data.copy() # Copy the variable data into memory
        #print("{} = {}".format(varname, var))
        
        # Look for coordinate variables for the plotting dimensions
        coord_vars = {}
        for dimname in dataset_var.dimensions:
            if (dimname in f.variables):
                coord_vars[dimname] = f.variables[dimname].data.copy()
        print("Coordinate variables found: {}".format(coord_vars))
        
        del dataset_var # Must be removed from scope for the dataset to close
        
    return var, var_dims, coord_vars

def nc_read(file, show_coordvars=False):
    """ Reads a simple format of NetCDF datasets with n dimensions (each having a 
    coordinate variable) and m (non-coordinate) variables. The dimension names and 
    lengths, the coordinate variables, the variable names and the variables are returned. 
    Note that variables are loaded into memory and returned as arrays, so this is not 
    adequate for very large datasets. """
    np.set_printoptions(linewidth=1000)
    
    with netcdf.netcdf_file(file, 'r') as f: # Open the NetCDF file
        print("Dimensions of the dataset: {}".format(f.dimensions)) # Print all dims
        
        # Define dictionaries to hold vars and their respective dimensions
        vars_dict = OrderedDict()
        dims_dict = OrderedDict()
        
        for el in f.variables:
            vars_dict[el] = f.variables[el].data.copy() # Copy the variable data into memory
            dims_dict[el] = f.variables[el].dimensions # Store the dimensions of each var
        
        # Look for coordinate variables for the plotting dimensions
        coord_vars = {}
        for dimname in f.dimensions:
            if (dimname in f.variables):
                coord_vars[dimname] = f.variables[dimname].data.copy()
        if (show_coordvars):
            print("Coordinate variables found: {}".format(coord_vars))
        else:
            print("Coordinate variables found: {}".format(coord_vars.keys()))
        
    return vars_dict, dims_dict, coord_vars

if __name__ == "__main__":
    if len(argv)>1:
        vars_dict, dims_dict, coord_vars = nc_read(argv[1], show_coordvars=False)
        print('\n{}'.format( vars_dict ))
        print('\n{}'.format( dims_dict ))
        print('\n{}'.format( coord_vars ))
