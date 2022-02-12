import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

from matplotlib import rc
from utopya import DataManager, UniverseGroup
from utopya.plotting import is_plot_func, PlotHelper, UniversePlotCreator

from .tools import titles

plt.rcParams.update({"text.usetex": True, "font.family": 'serif', "figure.figsize": [10, 8],
"font.size": 12, "font.weight": 'bold', "lines.linewidth": 3, "figure.dpi": 600})

# -----------------------------------------------------------------------------
def power_law(val, a, b, c: int = 0):
    return 1.0*a/((val+c)**b)

def Gauss(val, m, o, A):
    return A*np.exp(-1.0/o*(val-m)**2)

def var_ER(m, N):
    return m*(1-m/(N-1))

# -----------------------------------------------------------------------------
@is_plot_func(creator_type=UniversePlotCreator)
def degree_distribution(dm: DataManager, *,
                        uni: UniverseGroup,
                        hlpr: PlotHelper,
                        **plot_kwargs: dict):

    """Plots the degree distribution"""

    data = uni['data']['KronGen']['NetworkAnalyser']['graph_data']

    to_plot = 'degree_sequence' if 'degree_sequence' in data.keys() else 'degree'

    data_to_plot = np.asarray(data[to_plot].data)[0]

    if 'bins' in plot_kwargs.keys():
        plot_kwargs.remove('bins')

    # Plot histograms
    if (to_plot == 'degree_sequence'):
        x = np.nonzero(data_to_plot)[0]
        y = [data_to_plot[i] for i in x]
        num_vertices = np.sum(y)
        loc = np.sum(x*y)/num_vertices
        std = np.sqrt(np.sum(y*(x-loc)**2)/num_vertices)
        max = x[-1]
        min = x[0]
        hlpr.ax.scatter(x, y, **plot_kwargs)
    else:
        loc = np.mean(data_to_plot)
        std = np.std(data_to_plot)
        max = np.max(data_to_plot)
        min = np.min(data_to_plot)
        hist = hlpr.ax.hist(data_to_plot, **plot_kwargs)

    # Plot info box
    txt = (f"mean: {np.around(loc, 3)}"
        + r"$\pm$"+f"{np.around(std, 3)}"
        + "\n"
        + f"max: {np.around(max, 3)}"
        + "\n"
        + f"min: {np.around(min, 3)}")
    hlpr.ax.text(0.65, 0.85, txt, color="cornflowerblue",
                 transform=hlpr.ax.transAxes, backgroundcolor=(1, 1, 1, 0.8))
    hlpr.ax.axvline(loc, color="cornflowerblue", zorder=-1)
