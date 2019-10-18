# data_anal.py
""" Plots the specified variable as well as all other variables of the dataset defined on 
the same dimensions.
The following are specified as command-line arguments: (1) the path to the NetCDF 
dataset, (2) the name of the variable to be plotted, and (3) the comma-separated indices 
to be plotted for each of the variable's dimensions. (To specify which two dimensions to 
plot against, enter '-1' instead of an index.) """

import argparse
import numpy as np
import matplotlib.pyplot as plt

import plot_routines # Routine for plotting data
import nc_IO # Routine for reading data from a NetCDF file

# General functions
def indexing_tuple(input_dims, var_dims):
    """ Returns the tuple necessary for plotting on the desired dimensions and the names 
    of these dimensions. 'input_dims' is the plotting dimensions from the command line 
    and 'var_dims' is the list of dimensions on which the variable is defined. """
    
    # Check that the right number of dimensions were provided on the command line
    assert (len(var_dims) == len(input_dims)), "Incorrect number of dimensions was provided."
    
    # Build a tuple to index the variable
    plotting_dims = [] # This list will contain the names of the plotting dimensions
    for idx, val in enumerate(input_dims):
        if (val=='p'): # Want whole slice if value is negative
            input_dims[idx] = slice(None,None,None) # Equivalent to ":"
            plotting_dims.append(var_dims[idx])
        else:
            input_dims[idx] = int(input_dims[idx])
    tup = tuple(input_dims) # Redefine as tuple
    return tup, plotting_dims

def colorplor_extents(plotting_dims, coord_vars):
    assert (len(plotting_dims) == 2), "Function 'colorplor_extents' assumes two plotting dimensions."
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
    
    return horiz_extent, vert_extent

# Functions for the "multip" routines
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

# Plotting functions
def single_colorplot(plotting_dims, tup, coord_vars, varname, var):
    horiz_extent, vert_extent = colorplor_extents(plotting_dims, coord_vars)
    
    fig, ax = plt.subplots(1,1, subplot_kw={'aspect':'auto', 'adjustable':'box'})
    
    Labels = [plotting_dims[0], plotting_dims[1], varname] # Axis labels
    
    plot_routines.ColorPlot(fig, ax, var[tup], Labels, horiz_extent, vert_extent)
    ax.locator_params(axis='x', min_n_ticks=3) # Sets minimum tick number
    
    return

def single_lineplot(plotting_dims, tup, coord_vars, varname, var):
    # The horizontal axis for the plot
    if (plotting_dims[0] in coord_vars): # Get the coordinate variable if possible
        horiz_axis = coord_vars[plotting_dims[0]]
    else: # Otherwise just plot vs. index
        horiz_axis = range(len(var[tup]))
    
    fig, ax, cols = grid_plot(1,1, sharex='all')
    
    Labels = [plotting_dims[0], varname] # Axis labels
    
    plot_routines.Plot(ax, horiz_axis, var[tup], Labels)
    ax.locator_params(axis='x', min_n_ticks=3) # Sets minimum tick number
    
    return

def multip_colorplot(plotting_dims, tup, coord_vars, vars_to_plot, vars_dict, rows=2):
    horiz_extent, vert_extent = colorplor_extents(plotting_dims, coord_vars)
    
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
        
    return

def multip_lineplot(plotting_dims, tup, coord_vars, vars_to_plot, vars_dict, rows=2):
    # The horizontal axis for the plot
    if (plotting_dims[0] in coord_vars): # Get the coordinate variable if possible
        horiz_axis = coord_vars[plotting_dims[0]]
    else: # Otherwise just plot vs. index
        horiz_axis = range(len(var[tup]))
    
    fig, axes, cols = grid_plot(rows, len(vars_to_plot), sharex='all')
    for idx, val in enumerate(vars_to_plot):
        Labels = [plotting_dims[0], val] # Axis labels
        
        plot_routines.Plot(axes.flatten(order='F')[idx], horiz_axis, vars_dict[val][tup], Labels)
        axes.flatten(order='F')[idx].locator_params(axis='x', min_n_ticks=3) # Sets minimum tick number
    
    for row in range(rows): # Go through rows and cols to turn off x-axis labels except at bottom
        for col in range(cols):
            if (row!=rows-1):
                axes[row][col].set_xlabel('')
    return


def main(filename, varname, input_dims, same, display=True, save=True):
    print() # Linebreak
    
    # "coord_vars" is a dictionary containing coordinate variables
    vars_dict, dims_dict, coord_vars = nc_IO.nc_read(filename) # Get data
    var      = vars_dict[varname] # Get the requested variable
    var_dims = dims_dict[varname] # "var_dims" is a tuple with the dims on which the var is defined
    
    # Build a tuple to index the variable
    tup, plotting_dims = indexing_tuple(input_dims, var_dims)
    
    print("\nVariable dimensions: {}".format(var_dims))
    print(  "Plotting dimensions: {}".format(plotting_dims))
    # Check that either one or two dims were chosen for plotting
    assert (len(plotting_dims)==1 or len(plotting_dims)==2), "Number of plotting dimensions should be 1 or 2."
    
    if (same):
        # Get all the vars that share the same dimensions as "varname"
        vars_to_plot = find_vars(varname, dims_dict)
        
        if ( len(plotting_dims) == 2 ): # We plot every var in the dataset against the first two dimensions
            multip_colorplot(plotting_dims, tup, coord_vars, vars_to_plot, vars_dict)
        
        elif ( len(plotting_dims) == 1 ): # We plot every var in the dataset against the single plotting dimension
            multip_lineplot(plotting_dims, tup, coord_vars, vars_to_plot, vars_dict)
        
    else:
        if ( len(plotting_dims) == 2 ): 
            single_colorplot(plotting_dims, tup, coord_vars, varname, var)
        elif ( len(plotting_dims) == 1 ):
            single_lineplot(plotting_dims, tup, coord_vars, varname, var)
        
    plt.tight_layout()
    if (save):
        plt.savefig(filename+".pdf",bbox_inches='tight')
    if (display):
        plt.show()
    
    return


if __name__ == "__main__":
    # We parse the input using argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='Name of NetCDF dataset')
    parser.add_argument('variable', help='Name of variable to be plotted')
    parser.add_argument('dimensions', type=lambda s: [item for item in s.split(',')],
                        help='Comma-separated list of dimension indices to plot. Either one \
                        or two of these should be a "p", which indicates a plotting dimension.')
    parser.add_argument('-s', '--same', action='store_true', 
                        help='Plot all variables in the dataset defined on the same \
                        dimensions as "variable"')
    parser.add_argument('-m', '--mode', choices=['d','s','b'], default='d',
                        help='Set mode: "d" is for "display", "s" is for "save", and "b" is for "both".')
    
    args = parser.parse_args()
    
    # Determine whether the plot is to be displayed and/or saved
    display = True
    save    = True
    if (args.mode=='d'):
        save    = False
    elif (args.mode=='s'):
        display = False
    
    main(args.filename, args.variable, args.dimensions, args.same, display, save)
