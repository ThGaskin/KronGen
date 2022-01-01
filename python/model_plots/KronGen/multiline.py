import numpy as np
from utopya.plotting import is_plot_func, PlotHelper

# -----------------------------------------------------------------------------
@is_plot_func(use_dag=True, required_dag_tags=['data', 'the_dm'])
def multiline(data,
              hlpr: PlotHelper,
              x: str,
              target: str=None,
              plot_target_line: bool=True,
              **plot_kwargs):

    """For multiverse runs, this produces a line plot showing actual and any
    target values.

    Arguments:
        data (xarray): the dataset
        hlpr (PlotHelper): description
        x (str): the parameter dimension of the diagram.
        target (str): the parameter name
        plot_kwargs (dict, optional): kwargs passed to the pcolor plot function
    """

    to_plot = data['data']

    hlpr.ax.plot(to_plot.coords[x].data, to_plot.data_vars['y'].data, **plot_kwargs)

    hlpr.ax.grid(linewidth=0.5, alpha=0.5)

    # Find and plot any target values, if present
    dm = data['the_dm']._data['cfg']['run']._data['parameter_space']['KronGen']['create_graph']
    def finditem(key, obj):
        if key in obj: return obj[key]
        for k, v in obj.items():
            if isinstance(v,dict):
                item = finditem(key, v)
                if item is not None:
                    return item

    target = x if target is None else target
    target_vals = finditem(target, dm)
    if (target_vals is not None and plot_target_line):
        try:
            t = [i for i in target_vals]
        except:
            t = target_vals*np.ones(len(to_plot.coords[x].data))

        hlpr.ax.plot(to_plot.coords[x].data, t, color='black', linestyle='dashed',
                     label='Target', zorder=-1)
        hlpr.ax.legend()
