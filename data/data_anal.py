# data_anal.py
""" Plots the specified variable as well as all other variables of the dataset defined on 
the same dimensions.
The following are specified as command-line arguments: (1) the path to the NetCDF 
dataset, (2) the name of the variable to be plotted, and (3) the indices to be plotted 
for each of the variable's dimensions. (To specify which two dimensions to plot against, 
enter '-1' instead of an index.) """

from sys import argv # argv is the list of space-separated arguments from command prompt.
import numpy as np
import matplotlib.pyplot as plt

import plot_routines # Routine for plotting data
import nc_IO # Routine for reading data from a NetCDF file


def parse_args(args):
    """ Parses through the arguments to get the filename, varname, and dims. """
    # Unpack the argv list
    filename = argv[1] # Element 1 is dataset name (a string)
    varname  = argv[2] # Element 2 is variable name (a string)
    dims = argv[3:] # Subsequent elements are slices to plot. Build a list.
    for idx, val in enumerate(dims): # Make the list elements integers
        dims[idx] = int(val)
    
    return filename, varname, dims

def find_vars(varname, dims_dict):
    """ Finds all variables defined on the same dimensions as "var". """
    vars_to_plot = []
    for label in dims_dict:
        if (dims_dict[varname] == dims_dict[label]):
            vars_to_plot.append(label)
    
    return vars_to_plot

def grid_plot(rows, numplots, **kwargs):
    """ Outputs a grid of subplots given the number of subplots needed as well as the 
    desired number of rows. Also returns the number of columns. """
    # Plotting commands
    # Set the number of rows and columns in the plot
    if (numplots%rows==0):
        cols = numplots//rows
    else:
        cols = numplots//rows + 1 # Need extra column if remainder is nonzero
    
    # Create figure and subplots
    fig, axes = plt.subplots(rows, cols, subplot_kw={'aspect':'auto', 'adjustable':'box'}, figsize=(14.4,4.8), **kwargs)
    
    return fig, axes, cols


def main():
    filename, varname, dims = parse_args(argv)
    
    vars_dict, dims_dict, coord_vars = nc_IO.nc_read(filename) # Get data
    var      = vars_dict[varname] # Get the requested variable
    var_dims = dims_dict[varname] # Get its dimensions
    # "var_dims" is a tuple with the dims on which the var is defined
    # "coord_vars" is a dictionary containing coordinate variables
    
    # Check that the right number of dimensions were provided on the command line
    assert (len(var_dims) == len(dims)), "Incorrect number of dimensions was provided."
    
    
    # Build a tuple to index the variable
    tup = dims # Starts off as a list
    plotting_dims = [] # This list will contain the names of the plotting dimensions
    for idx, val in enumerate(tup):
        if (val<0): # Want whole slice if value is negative
            tup[idx] = slice(None,None,None) # Equivalent to ":"
            plotting_dims.append(var_dims[idx])
    tup = tuple(tup) # Redefine as tuple
    
    
    print("Plotting dimensions: {}".format(plotting_dims))
    # Check that either one or two dims were chosen for plotting
    assert (len(plotting_dims)==1 or len(plotting_dims)==2), "Number of plotting dimensions should be 1 or 2."
    
    # Get all the vars that share the same dimensions as "varname"
    vars_to_plot = find_vars(varname, dims_dict)
    
    if ( len(plotting_dims) == 2 ): # We plot every var in the dataset against the first two dimensions
        
        # Define the extents of the axes if they have coordinate variables
        if (plotting_dims[0] in coord_vars):
            X = coord_vars[plotting_dims[0]]
            dx = (X[1] - X[0])/2.
            horiz_extent = (X[0]-dx, X[-1]+dx)
        else:
            horiz_extent = () # Leads to default behaviour
        
        if (plotting_dims[1] in coord_vars):
            Y = coord_vars[plotting_dims[1]]
            dy = (Y[1] - Y[0])/2.
            vert_extent = (Y[0]-dy, Y[-1]+dy)
        else:
            vert_extent = () # Leads to default behaviour
        
        
        rows = 2 # Number of rows of plots
        fig, axes, cols = grid_plot(rows, len(vars_to_plot), sharex='all', sharey='all')
        for idx, val in enumerate(vars_to_plot):
            Labels = [plotting_dims[0], plotting_dims[1], val] # Axis labels
            plot_routines.ColorPlot(fig, axes.flatten(order='F')[idx], vars_dict[val][tup], Labels, horiz_extent, vert_extent)
            axes.flatten(order='F')[idx].locator_params(axis='x', min_n_ticks=3) # Sets minimum tick number
        
        for row in range(rows): # Go through rows and cols to turn off axis labels except at edges
            for col in range(cols):
                if (row!=rows-1):
                    axes[row][col].set_xlabel('')
                if (col!=0):
                    axes[row][col].set_ylabel('')
        
        plt.tight_layout()
        plt.savefig(filename+".pdf",bbox_inches='tight')
        plt.show()
    
    elif ( len(plotting_dims) == 1 ): # We plot every var in the dataset against the single plotting dimension
        
        # The horizontal axis for the plot
        if (plotting_dims[0] in coord_vars): # Get the coordinate variable if possible
            horiz_axis = coord_vars[plotting_dims[0]]
        else: # Otherwise just plot vs. index
            horiz_axis = range(len(var[tup]))
        
        
        rows = 2 # Number of rows of plots
        fig, axes, cols = grid_plot(rows, len(vars_to_plot), sharex='all')
        for idx, val in enumerate(vars_to_plot):
            Labels = [plotting_dims[0], val] # Axis labels
            #plot_routines.ColorPlot(fig, axes.flatten(order='F')[idx], vars_dict[val][tup], Labels, horiz_extent, vert_extent)
            plot_routines.Plot(axes.flatten(order='F')[idx], horiz_axis, vars_dict[val][tup], Labels)
            axes.flatten(order='F')[idx].locator_params(axis='x', min_n_ticks=3) # Sets minimum tick number
        
        for row in range(rows): # Go through rows and cols to turn off x-axis labels except at bottom
            for col in range(cols):
                if (row!=rows-1):
                    axes[row][col].set_xlabel('')
        
        plt.tight_layout()
        plt.savefig(filename+".pdf",bbox_inches='tight')
        plt.show()
    
    return

if __name__ == "__main__":
    main()
