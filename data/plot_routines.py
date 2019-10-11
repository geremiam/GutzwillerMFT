# plot_routines.py
""" Various routines for plotting specific data formats. """

import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import numpy as np

def ColorPlot(fig, ax, data, Labels, horiz_extent=(), vert_extent=(), logscale=False):
    """ Plots the 2D array "data" (first index along horizontal axis and second index 
    along vertical axis). The three strings in the list "Labels" correspond to labels for 
    the horizontal, vertical, and color axes, respectively. Axis ends can be specified as 
    tuples (optional). """
    
    # Check for axis extents
    if (horiz_extent==()):
        horiz_extent = (-0.5, data.shape[0]-0.5) # Default behaviour
    if (vert_extent==()):
        vert_extent = (-0.5, data.shape[1]-0.5) # Default behaviour
    extent = horiz_extent + vert_extent
    
    
    cmaps = ['PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 
             'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']
    cmap=cmaps[2]
    
    # Center the colormap at zero.
    vmax = np.amax(np.abs(data))
    vmin = -vmax
    
    if (logscale):
        norm = SymLogNorm(1.e-6,vmin=vmin,vmax=vmax)
    else:
        norm = None
    
    # Transpose the data because imshow uses 'ij' indexing instead of 'xy'
    data = np.transpose(data)
    im = ax.imshow(data, cmap=cmap, aspect='auto', interpolation='nearest', 
                   vmin=vmin, vmax=vmax,origin='lower', extent=extent, norm=norm)
    
    cbar = fig.colorbar(im, ax=ax, format='%g', pad=0.03)
    
    #ax.contour(tensor, 15, colors='k', vmin=vmin, vmax=vmax,origin='lower', extent=extent)
    #ax.contour(tensor, (0,), colors='C3', vmin=vmin, vmax=vmax, origin='lower', extent=extent, linewidths=2)
    
    # Set axis labels
    ax.set_xlabel(Labels[0])
    ax.set_ylabel(Labels[1])
    cbar.set_label(Labels[2])
    
    return

def Plot(ax, x, y, Labels):
    """ Plots the 2D array "data" (first index along horizontal axis and second index 
    along vertical axis). The three strings in the list "Labels" correspond to labels for 
    the horizontal, vertical, and color axes, respectively. Axis ends can be specified as 
    tuples (optional). """
    
    ax.plot(x, y, color='tab:grey', zorder=1)
    ax.scatter(x, y, color='tab:red', marker='.', s=16., zorder=2)
    
    # Set axis labels
    ax.set_xlabel(Labels[0])
    ax.set_ylabel(Labels[1])
    
    return
