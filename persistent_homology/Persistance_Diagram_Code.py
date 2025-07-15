#importing packages
import sys
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

#Modifications to the plotting
from os import path
from math import isfinite
import numpy as np
from functools import lru_cache

__author__ = "Vincent Rouvreau, Bertrand Michel, Theo Lacombe"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

_gudhi_matplotlib_use_tex = True

##setting up the visuals of the plot
def __min_birth_max_death(persistence, band=0.0):
    """This function returns (min_birth, max_death) from the persistence.

    :param persistence: The persistence to plot.
    :type persistence: list of tuples(dimension, tuple(birth, death)).
    :param band: band #!!What does band mean? What is it keeping track of?
    :type band: float.
    :returns: (float, float) -- (min_birth, max_death).
    """
    # Look for minimum birth date and maximum death date for plot optimisation
    max_death = 0
    min_birth = persistence[0][1][0] #!!What do the numbers in brackets mean?
    for interval in reversed(persistence):
        if float(interval[1][1]) != float("inf"):
            if float(interval[1][1]) > max_death:
                max_death = float(interval[1][1])
        if float(interval[1][0]) > max_death:
            max_death = float(interval[1][0])
        if float(interval[1][0]) < min_birth:
            min_birth = float(interval[1][0]) #Why do they have the two for max_birth but only the one for min_birth
    if band > 0.0:
        max_death += band #what does the += mean
    return (min_birth, max_death)

def _array_handler(a):
    '''
    :param a: if array, assumes it is a (n x 2) np.array and return a
                persistence-compatible list (padding with 0), so that the
                plot can be performed seamlessly.
    '''
    if isinstance(a[0][1], np.float64) or isinstance(a[0][1], float):
        return [[0, x] for x in a]
    else:
        return a
    
@lru_cache(maxsize=1)

def _matplotlib_can_use_tex():
    """This function returns True if matplotlib can deal with LaTeX, False otherwise.
    The returned value is cached.
    """
    try:
        from matplotlib import checkdep_usetex
        return checkdep_usetex(True)
    except ImportError:
        print("This function is not available, you may be missing matplotlib.")

def plot_persistence_diagram(
    persistence=[],
    alpha=0.6, #transparency of the graph
    band=0.0,
    max_intervals=1000, #max view 1000 intervals with the longest lifetime
    inf_delta=0.1, #what does this do?
    legend=False,
    colormap=None,
    axes=None,
    fontsize=16,
    greyblock=True,
    use_fixed_max = True
) :

    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib import rc
        if _gudhi_matplotlib_use_tex and _matplotlib_can_use_tex():
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif')
        else:
            plt.rc('text', usetex=False)
            plt.rc('font', family='DejaVu Sans')

        persistence = _array_handler(persistence)

        if max_intervals > 0 and max_intervals < len(persistence):
            # Sort by life time, then takes only the max_intervals elements
            persistence = sorted(
                persistence,
                key=lambda life_time: life_time[1][1] - life_time[1][0],
                reverse=True,
            )[:max_intervals]

        if colormap == None:
            colormap = plt.cm.Set1.colors
        if axes == None:
            fig, axes = plt.subplots(1, 1)
        
        (min_birth, max_death) = __min_birth_max_death(persistence, band)
        delta = (max_death - min_birth) * inf_delta
        # Replace infinity values with max_death + delta for diagram to be more
        # readable
        if use_fixed_max: #!!why does it create points with a fixed max? And where is it getting the 'use_fixed_max' from
            infinity = 7942/60 + delta
            axis_end = 7942/60 + delta / 2
        else:
            infinity = max_death + delta
            axis_end = max_death + delta / 2
#         axis_start = min_birth - delta
        axis_start = 0

        # bootstrap band
        if band > 0.0:
            x = np.linspace(axis_start, infinity, 1000)
            axes.fill_between(x, x, x + band, alpha=alpha, facecolor="red")
        # lower diag patch
        if greyblock:
            axes.add_patch(mpatches.Polygon([[axis_start, axis_start], [axis_end, axis_start], [axis_end, axis_end]], fill=True, color='lightgrey'))
        # Draw points in loop
        pts_at_infty = False  # Records presence of pts at infty #!!How is a point at infinity? What is considered infinity?
        for interval in reversed(persistence):
            if float(interval[1][1]) != float("inf"):
                # Finite death case
                if interval[0]==0:
                    axes.scatter(
                        interval[1][0],
                        interval[1][1],
                        alpha=alpha,
                        color=colormap[interval[0]],
                    )
                else:
                    axes.scatter(
                        interval[1][0],
                        interval[1][1],
                        alpha=alpha,
                        color=colormap[interval[0]],
                        marker = 'X'
                    )
            else:
                pts_at_infty = True
                # Infinite death case for diagram to be nicer. !!What is an infinite death case
                if interval[0]==0:
                    axes.scatter(
                        interval[1][0], infinity, alpha=alpha, color=colormap[interval[0]]
                    )
                else:
                    axes.scatter(
                        interval[1][0], infinity, alpha=alpha, color=colormap[interval[0]], marker = 'X'
                    )
        if pts_at_infty: #Creating the inifinity line at the top of the persistance diagram
            axes.plot([axis_start, axis_end], [axis_start, axis_end], linewidth=1.0, color="k")
            axes.plot([axis_start, axis_end], [infinity, infinity], linewidth=1.0, color="k", alpha=alpha)
            # Infinity label
            yt = axes.get_yticks()
            yt = yt[np.where(yt < axis_end)] # to avoid ploting ticklabel higher than infinity
            yt = np.append(yt, infinity)
            ytl = ["%.0f" % e for e in yt]  # to avoid float precision error
            ytl[-1] = r'$+\infty$'
            axes.set_yticks(yt)
            axes.set_yticklabels(ytl)

        if legend:
            dimensions = list(set(item[0] for item in persistence))
            axes.legend(
                handles=[
                    mpatches.Patch(color=colormap[dim], label=str(dim) + "D homology class")
                    for dim in dimensions
                ]
            )  

            axes.set_xlabel("Birth", fontsize=fontsize)
            axes.set_ylabel("Death", fontsize=fontsize)
            axes.set_title("Persistence diagram", fontsize=fontsize)
            # Ends plot on infinity value and starts a little bit before min_birth
            axes.axis([axis_start, axis_end, axis_start, infinity + delta/2])
            return axes
        
    except ImportError:
        print("This function is not available, you may be missing matplotlib.")
        

#importing stuff again to set up data and simplex formation
import numpy as np
import gudhi as gd
import gudhi.weighted_rips_complex #see if we can exclude this one

def get_tda_info(citypath, city, max_hom_dim = 1):
    #are we supposed to upload/have our own citypath, as this is not defined later in the code. Same with city, do we create a vector with them?
    D_city = np.load(citypath + f'/{city}_d_matrix.npy')
    #D_city = np.genfromtxt(citypath + f'/{city}_d_matrix.csv', delimiter = ",")
    wait_city = np.genfromtxt(citypath + f'/{city}_waits.csv', delimiter = ",")
     #So then, is the citypath referring to the distance matrix for each of the cities? Also do not have to upload the waittimes.
        #May not have to copy this specific code but just relocate files to folder with code

    # Calculate simplex pairs for the homology classes
    cpx = gd.weighted_rips_complex.WeightedRipsComplex(distance_matrix = D_city, weights = wait_city).create_simplex_tree(max_dimension = max_hom_dim + 1)
        #Have to figure out how to alter this code as not to get weighted complex
    cpx.compute_persistence()
    ph = cpx.persistence()
    all_pairs = cpx.persistence_pairs()
    return cpx, ph, all_pairs

#plotting the persistance diagrams
import geopandas as gpd
from shapely.geometry import LineString, Polygon
import scipy.stats

def plot_deathsimplices(citypath, city, cpx, all_pairs, hom_dim, ax, criteria = 'death', zscore_thresh = 1, citydf_fname = None, vmin0 = 3235, vmin1 = 4410, vmax0 = 7380, vmax1 = 7942, legend = False):
#Where do they get the end values, including the zscore threshold of 1 (why not 1.05). Or, how did they find the highest death filtration of the 0th 1st homology

        # Load geographic data
    if citydf_fname is None:
        city_df = gpd.read_file(citypath + f'/{city}_zip.geojson')
    else:
        city_df = gpd.read_file(citypath + "/" + citydf_fname)
    polls_df = gpd.read_file(citypath + f'/{city}_polls.geojson')     
