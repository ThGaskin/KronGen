import matplotlib.pyplot as plt
import numpy as np
import math

from utopya import DataManager, UniverseGroup
from utopya.plotting import is_plot_func, PlotHelper, UniversePlotCreator

# -----------------------------------------------------------------------------
@is_plot_func(creator_type=UniversePlotCreator)
def network_stats(dm: DataManager, *,
                  uni: UniverseGroup,
                  hlpr: PlotHelper,
                  **plot_kwargs: dict):

    """Plots a sheet with the network topology stats"""

    data = uni['data']['KronGen']['NetworkAnalyser']['graph_data']
    plots = list(data.keys())[2:] # discard vertices and edges
    n_plots = len(plots)

    # .. Setup the figure with two columns and as many rows as necessary .......
    n_rows = math.ceil(n_plots/2)
    height_ratios = np.ones(n_rows)
    figure = plt.figure(figsize=(10, n_rows * 5))
    axs = []
    gs = figure.add_gridspec(ncols=2, nrows=n_rows,
                  height_ratios=height_ratios, width_ratios=[1, 1], hspace=0.2)
    gridspec = []
    plots_added = 0
    for i in range(n_plots):
        for j in range(2):
            if (plots_added < n_plots):
                gridspec.append((i, j))
                plots_added+=1
            else:
                break
    for item in gridspec:
        axs.append([figure.add_subplot(gs[item])])
    hlpr.attach_figure_and_axes(fig=figure, axes=axs)

    # .. Plot data histogram on each axis ......................................
    for i in range(n_plots):
        hlpr.select_axis(0, i)
        data_to_plot = np.asarray(data[plots[i]].data)[0]
        hist = hlpr.ax.hist(data_to_plot, **plot_kwargs)
        if (plots[i] == 'distance_max'):
            loc = np.max(data_to_plot)
            txt = f"max: {np.around(loc, 3)}"
            hlpr.ax.text(0.89*loc, np.max(hist[0]), txt, color="red")
        else:
            loc = np.mean(data_to_plot)
            txt = f"mean: {np.around(loc, 3)}"
            hlpr.ax.text(1.01*loc, np.max(hist[0]), txt, color="red")
        hlpr.ax.axvline(loc, color="red")
        hlpr.ax.set_title(plots[i])
        hlpr.ax.set_yscale('log')
