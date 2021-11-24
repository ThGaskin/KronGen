import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import csv
import pandas as pd
from matplotlib.lines import Line2D

from utopya.plotting import is_plot_func, PlotHelper

# -----------------------------------------------------------------------------
@is_plot_func(use_dag=True, required_dag_tags=['data'])
def multiline(data,
              hlpr: PlotHelper,
              x: str = None,
              y: list = None,
              **plot_kwargs):

    hlpr.ax.grid(linewidth=0.5, alpha=0.5)

    to_plot = data['data']

    # hlpr.ax.plot(to_plot.coords[x].data, to_plot.data)
    #
    #
    # df = pd.DataFrame(to_plot)
    # df.to_csv("~/test.csv")

    # Get x coordinates
    if x is None:
        for key in to_plot.coords.keys():
            if (key != 'vertex_idx'):
                x = key

    y_data = y if y is not None else to_plot.data_vars.keys()

    # Extract additional kwargs for specific lines
    plot_kwargs_separate = {}
    for var in y_data:
        if var in plot_kwargs:
            plot_kwargs_separate[var] = plot_kwargs[var]
            plot_kwargs.pop(var)
    
    # Get y coordinates and plot
    for var in y_data:
        additional_kwargs = {}
        if var in plot_kwargs_separate.keys():
            additional_kwargs.update(plot_kwargs_separate[var])

        hlpr.ax.plot(to_plot.coords[x].data, to_plot.data_vars[var].data, **plot_kwargs, **additional_kwargs)

    hlpr.ax.legend([Line2D([0], [0], color='black', linestyle='dashed', lw=1)],
                   ['Target value'])
