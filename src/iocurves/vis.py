#coding=utf-8
"""
Visualisations for the input-output function of neurons over time.
"""
import logging
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cm, colors, pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from src.config import settings
from src.io_time_heatmap import fig_time, save_location
from src.iocurves.analysis import get_data, get_inst_firing_rate, get_params, get_var

logger = logging.getLogger("io curves vis")


def plot_compare_dcl(title, results, ifr_windowsize=0.01, combine=False):
    """
    Plot voltage, chloride, and number spikes over time.

    :param title: Figure title
    :type title: str
    :param results:
    :type results: list of pd.DataFrame
    :param ifr_windowsize: Instantanenous window size (in seconds)
    :type ifr_windowsize: float
    :param combine: Display static and dynamic cl- traces on a single axis (True) or on their own axes (False)
    :type combine: bool
    :return: Figure and axes create in method
    :rtype: (plt.Figure, plt.Axes)
    """
    f, axes = plt.subplots(3, 1 + (1 - int(combine)), sharex='all', sharey='row')
    f.suptitle(title)
    prev_legend = []
    for i, result in enumerate(results):
        multiIndex = result.columns
        recorded_vars = multiIndex.levels[1]
        for i_var, recorded_var in enumerate(recorded_vars):
            legend = False
            ax = axes[i_var] if combine else axes[i_var, i]
            style = '--' if combine and i == 0 else '-'

            trial_results, mean = get_var(result, recorded_var)
            if recorded_var == 'spikes':
                # get number of spikes in total
                legend = [mean.iloc[-1]]
                if combine:
                    name = 'static' if i == 0 else 'dynamic'
                    prev_legend.append("{} ({})".format(legend[0], name))
                # convert spike times to instantaneous firing rate
                trial_length = trial_results.shape[1] + 1
                for j in range(1, trial_length):
                    trial_results.loc[:, j] = get_inst_firing_rate(trial_results[j], time_bin=ifr_windowsize)
                mean.loc[:] = np.mean(trial_results, axis=1)

            trial_results.plot(ax=ax, color='k', alpha=0.1, style=style, legend=legend)
            mean.plot(ax=ax, color='k', alpha=0.8, style=style, legend=legend)
            if legend and combine:
                ax.legend(prev_legend, loc='lower right')
            elif legend:
                ax.legend(legend, loc='lower right')

            if i == 0:
                # first column (Y axis labels)
                if recorded_var == 'v':
                    ylabel = 'Membrane Potential (mV)'
                elif recorded_var == 'cli':
                    ylabel = '$[Cl^{-}]_i$ (mM)'
                elif recorded_var == 'ifr' or recorded_var == 'spikes':
                    ylabel = "Instantaneous Firing Rate (Hz) \n "\
                             "[window size of {} ms]".format(ifr_windowsize*1000)
                else:
                    ylabel = None
                ax.set_ylabel(ylabel)

    if combine:
        axes[0].set_title("Static Chloride vs Dynamic Chloride")
    else:
        axes[0, 0].set_title("Static Chloride")
        axes[0, 1].set_title("Dynamic Chloride")

    for axis in axes.flatten():
        sns.despine(ax=axis)

    plt.xlabel('Time (ms)')
    # plt.xlim([0, (len(t_rec_) - 1) * 0.025])
    return f, axes


def io_curve(io_runs, ifr_windowsize=0.1, combine=False, heatmap=False, fill=False, time_points=None,
             show_cli=False, save_args=None):
    """
    Plot input-output of a neuron.

    Can be:

    * heatmap and combine - plot static and dynamic heatmaps AND the difference between them
    * heatmap - plot static and dynamic heatmaps
    * combine - plot static and dynamic input-output curves on one axes (like `plot_compare_dcl`)
    * neither - plot static and dynamic input-output curves on separate axes (like `plot_compare_dcl`)
    * io_runs is a dictionary - plot a specific run. Only if not heatmap and not combine.

    :param io_runs: Static and dynamic chloride dataframes with firing rate info (from `sim`).
    :type io_runs: list of pd.DataFrame or dict of pd.DataFrame
    :param ifr_windowsize: Window to compute instantanenous firing rate (in seconds) .
    :type ifr_windowsize: float
    :param combine: Plot static and dynamic on one axes or create a difference heatmap.
    :type combine: bool
    :param heatmap: Plot as heatmap(s).
    :type heatmap: bool
    :param fill: (Forward) Fill the dataframe.
    :type fill: bool
    :param time_points: Time points to plot (or None for a 3D scatter plot)
    :type time_points: list of float
    :param show_cli: Include internal chloride ion concentration heatmap (when heatmap is True, ignored otherwise).
    :type show_cli: bool
    :param save_args: File name (`fig_name`) and file types (`formats`) for saving.
    :type save_args: dict
    :return: Pair of figure, ax used
    :rtype: (plt.Figure, plt.Axes)
    """

    subplot_kw = None
    gridspec_kw = None
    if time_points is None:
        subplot_columns = 1
        heatmap = True
        subplot_kw = {"projection": '3d'}
    else:
        subplot_columns = len(time_points) + 1
        gridspec_kw = {'width_ratios': [5]*len(time_points) + [1]}

    extra = ''
    y_label_cl = None
    heatmap_kwargs = {
        "square": True,
        "vmin":   0,
        "vmax":   120,
        }
    label = 'IFR (Hz)'
    if heatmap and combine:
        f, axes = plt.subplots(3 + show_cli, subplot_columns, figsize=(settings.PAGE_W_FULL,
                                                                       settings.PAGE_H_half + settings.PAGE_H_4th),
                               sharex='col', sharey=False,
                               subplot_kw=subplot_kw, gridspec_kw=gridspec_kw)
        dfs = []
        index = None
        exc, inh = None, None
        ylabel = 'Relative\ninhibitory\nconductance'
        for i_cl, cl_state in enumerate(io_runs):
            params = get_params(list(cl_state.keys())[0])
            file_name = params[0]
            line = '-' if 'KCC2' in file_name else '--'
            axes_subplot = axes[i_cl] if time_points is None else axes[i_cl, :]
            FRdf, exc, inh = get_data(cl_state, ifr_windowsize, time_points)
            _plot_data(FRdf, exc, inh, label, axes_subplot, f, heatmap, fill, time_points, line,
                       first_plot=(i_cl == 0), cbar_ax=axes_subplot[-1], heatmap_kwargs=heatmap_kwargs)
            dfs.append(FRdf)
            y_label_cl = settings.DYNAMIC_CHLORIDE_STR_ABBR if 'KCC2' in file_name else\
                settings.STATIC_CHLORIDE_STR_ABBR
            axes_subplot = axes[i_cl] if time_points is None else axes[i_cl, 0]
            axes_subplot.set_ylabel(y_label_cl + "\n" + ylabel, rotation=0, ha='right', va='center_baseline')
            if show_cli and 'KCC2' in file_name:
                if show_cli:
                    axes_subplot = axes[-1] if time_points is None else axes[-1, :]
                    # heatmap_kwargs["vmax"] = vmax
                    cli_df, exc, inh = get_data(cl_state, ifr_windowsize, time_points, var='cli')
                    cmap = sns.cubehelix_palette(16, start=2.7, rot=-.2, light=0.98, dark=0.40, as_cmap=True)
                    from matplotlib.colors import LogNorm
                    cl_min = cli_df.values.min()
                    cl_max = cli_df.values.max()
                    _plot_data(cli_df, exc, inh, f"{settings.CLI} (mM)", axes_subplot, f, heatmap, fill,
                               time_points, line, first_plot=False, cmap=cmap, cbar_ax=axes_subplot[-1],
                               cbar_kws={'ticks': np.logspace(np.log10(cl_min), np.log10(cl_max), 6, base=10)},
                               heatmap_kwargs={
                                   'square': heatmap_kwargs['square'],
                                   "vmin":   cl_min,
                                   "vmax":   cl_max,
                                   "norm":   LogNorm(vmin=cl_min, vmax=cl_max)
                                   })
        logger.info("plotting difference")
        # do extra 'Difference' plot
        dif_df = dfs[1] - dfs[0]
        # get maximum value of dataframe
        vmax = dif_df.values.max()
        axes_subplot = axes[2] if time_points is None else axes[2, :]
        heatmap_kwargs["vmax"] = vmax
        _plot_data(dif_df, exc, inh, f"{settings.DELTA}{label}", axes_subplot, f, heatmap, fill, time_points,
                   first_plot=False,
                   cmap='viridis', cbar_ax=axes_subplot[-1],
                   heatmap_kwargs=heatmap_kwargs)
        y_label_cl = "Difference"
        axes_subplot = axes[2] if time_points is None else axes[2, 0]
        axes_subplot.set_ylabel(y_label_cl + "\n" + ylabel, rotation=0, ha='right', va='center_baseline')
        if time_points is None:
            ax = axes[-1]
        else:
            ax = axes[-1, 0]
        f.subplots_adjust(wspace=0.1, hspace=0.2, left=0.3, right=0.93)
    elif type(io_runs) == dict:
        # just plot the one
        f, axes = plt.subplots(1, subplot_columns, sharey='row', subplot_kw=subplot_kw)
        params = get_params(list(io_runs.keys())[0])
        file_name = params[0]
        line = '-' if 'KCC2' in file_name else '--'
        FRdf, exc, inh = get_data(io_runs, ifr_windowsize, time_points)
        _plot_data(FRdf, exc, inh, label, axes, f, heatmap, fill, time_points, heatmap_kwargs=heatmap_kwargs)

        ax = axes[0]
        extra = ' Dynamic' if 'KCC2' in file_name else ' Static'
    elif combine:
        # plot both cl_states on one set of axes
        f, axes = plt.subplots(1, subplot_columns, sharey='row', subplot_kw=subplot_kw)
        for i_cl, cl_state in enumerate(io_runs):
            params = get_params(cl_state.keys()[0])
            file_name = params[0]
            line = '-' if 'KCC2' in file_name else '--'
            FRdf, exc, inh = get_data(cl_state, ifr_windowsize, time_points)
            _plot_data(FRdf, exc, inh, label, axes, f, heatmap, fill, time_points, line,
                       first_plot=(i_cl == 0), cbar_ax=axes[-1],
                       heatmap_kwargs=heatmap_kwargs)
        ax = axes[0]
    else:
        # plot cl_states on multiple axes
        f, axes = plt.subplots(2, subplot_columns, sharex='all', sharey='all', subplot_kw=subplot_kw)
        for i_cl, cl_state in enumerate(io_runs):
            params = get_params(cl_state.keys()[0])
            file_name = params[0]
            line = '-' if 'KCC2' in file_name else '--'
            FRdf, exc, inh = get_data(cl_state, ifr_windowsize, time_points)
            _plot_data(FRdf, exc, inh, label, axes[i_cl, :], f, heatmap, fill, time_points, line,
                       first_plot=(i_cl == 0), cbar_ax=axes[-1],
                       heatmap_kwargs=heatmap_kwargs)
            if heatmap:
                y_label_cl = 'Dynamic ' if 'KCC2' in file_name else 'Static '
                axes[i_cl, 0].set_ylabel(y_label_cl + 'Inhibition')
        ax = axes[1, 0]

    title = list(io_runs.keys())[0] if type(io_runs) == dict else list(io_runs[0].keys())[0]
    f.suptitle(title + "\n" + "Input-Output curve" + extra)

    # fit to a sigmoid curve
    # xdata = np.log10(FRdf.index[1:].values)
    # ydata = FRdf.iloc[1:,0].values
    # popt, pcov = curve_fit(sigmoid, xdata, ydata)
    # x = np.linspace(0, np.log10(FRdf.index[-1]), 100)
    # y = sigmoid(x, *popt)
    # ax[i].semilogx(10**xdata,ydata,'o',10**x,y,'-')

    # ax[i].set_xlabel(str(time_point) + 'ms')
    if heatmap and y_label_cl is None:
        ax.set_ylabel("Relative inhibitory conductance")
    elif not heatmap:
        ax.set_ylabel("IFR (Hz)")
        ax.legend(title='Inhibition')
    ax.set_xlabel("Relative excitatory conductance")

    if save_args is not None and type(save_args) is dict:
        for fig_format in save_args["formats"]:
            f.savefig("{save_location}{fname}_{fig_time}_{window_size}{combine}{heatmap}{fill}.{format}"
                      .format(save_location=save_location,
                              fname=save_args["fig_name"],
                              fig_time=fig_time,
                              window_size=ifr_windowsize*1000,
                              combine="_combine" if combine else "",
                              heatmap="_heatmap" if heatmap else "",
                              fill="_fill" if fill else "",
                              format=fig_format))

    return f, axes


def _plot_data(FRdf, exc, inh, label, ax, fig, heatmap=False, fill=False, time_points=None, line='-',
               first_plot=True, cmap=None, cbar_ax=None, cbar_kws=None, heatmap_kwargs=None):
    """
    Plot data from the dataframe as heatmap, input-output curves, or 3D scatter plot (not advised).

    :param FRdf: Firing rate data.
    :type FRdf: pd.DataFrame
    :param exc: List of excitation values in FRdf.
    :type exc: list
    :param inh: List of inhibition values in FRdf.
    :type inh: list
    :param label: Display name for FRdf values (e.g. Firing Rate or IFR).
    :type label: str
    :param ax: Axes to plot on. If t_points is not None, must be the same-size array.
    :type ax: np.ndarray[plt.Axes] or plt.Axes or Axes3D
    :param fig: Figure to used to plot colorbar if time_pints is None.
    :type fig: plt.Figure
    :param heatmap: Display i-o curves as heatmap instead of multiple lines.
    :type heatmap: bool
    :param fill: Fill empty values of FRdf ('ffill' used).
    :type fill: bool
    :param time_points: Time points to take from FRdf to plot.
    :type time_points: list
    :param line: Style of line if heatmap is False.
    :type line: str
    :param first_plot: Add the time point as a title. Ignored if time_points is None.
    :type first_plot: bool
    :param cmap: Colormap to use (defaults to matplotlib's default of 'viridis')
    :type cmap: str
    :param cbar_ax: Where to plot the colorbar when time_points is not None. If None, default behaviour is to take
        space from axes for each time point.
    :type cbar_ax: plt.Axes
    :param cbar_kws: Additional keywords for colorbar.
    :type cbar_kws: dict
    :param heatmap_kwargs: Additional keywords for `sns.heatmap`
    :type heatmap_kwargs: dict
    """

    # ensure there is a MultiIndex to reference
    index = pd.MultiIndex.from_product([inh, exc], names=['in', 'ex'])
    if cbar_kws is None:
        cbar_kws = dict()

    if time_points is None:
        _plot_3d_scatter(FRdf, ax, cmap, fig, fill, heatmap_kwargs, index, label)
    else:
        for i, time_point in enumerate(time_points):
            df_time_point = FRdf.loc[time_point].reindex(index).sort_index()
            if heatmap:
                _plot_io_heatmap(ax[i], cbar_ax, cbar_kws, cmap, fill, heatmap_kwargs, label, df_time_point,
                                 i == 0, i == len(time_points) - 1)
            else:
                if fill:
                    df_time_point = df_time_point.fillna(method='ffill')
                for inh in df_time_point.index.levels[0]:
                    df_time_point.loc[inh].plot(ax=ax[i], logx=True, linestyle=line, label=str(inh))

            if first_plot:
                ax[i].set_title(str(time_point) + 'ms')


def _plot_io_heatmap(ax, cbar_ax, cbar_kws, cmap, fill, heatmap_kwargs, label, df_time_point,
                     keep_tick_labels, add_colorbar):
    """
    Plot heatmap of excitation & inhibition vs firing rate (in ti).

    :param ax: Axes for plotting
    :type ax: plt.Axes
    :param cbar_ax: Axes for colorbar (if idx
    :type cbar_ax: plt.Axes
    :param cbar_kws: Keywords for colorbar (except label)
    :type cbar_kws: dict
    :param cmap: Color to use for colorbar
    :type cmap: str or Colormap
    :param fill: (Forward) Fill in nan values of heatmap.
    :type fill: bool
    :param heatmap_kwargs: Keywords for `sns.heatmap`
    :type heatmap_kwargs: dict
    :param label: Colorbar label
    :type label: str
    :param df_time_point:
    :type df_time_point: pd.DataFrame
    :param keep_tick_labels: First column has ticklabels
    :type keep_tick_labels: bool
    :type add_colorbar: Last column creates colorbar in cbar_ax (if it's not None)
    :type add_colorbar: bool
    """
    heatmap_df = pd.DataFrame()
    # reformat Dataframe so it's easier to fill and plot
    for (inh, exc) in df_time_point.index:
        heatmap_df.loc[exc, inh] = df_time_point.loc[(inh, exc)]
    # have exc on X axis (.T) and reverse inh (iloc[::-1])
    heatmap_df = heatmap_df.T.iloc[::-1]
    if fill:
        heatmap_df = heatmap_df.fillna(method='ffill')
    # create colorbar for every heatmap if cbar_ax is None, otherwise only add colorbar for last coloumn
    cbar = True if cbar_ax is None else add_colorbar
    if "ticks" in cbar_kws and cbar_ax is not None and cbar:
        ticks = [int(t) for t in cbar_kws["ticks"]]
    hm = sns.heatmap(heatmap_df, ax=ax,
                     cbar=cbar, cbar_ax=cbar_ax,
                     cbar_kws=dict(label=label, **cbar_kws),
                     cmap=cmap,
                     **heatmap_kwargs)
    if "ticks" in cbar_kws and cbar_ax is not None and cbar:
        cbar_ax.yaxis.set_ticklabels(ticks)
        cbar_ax.yaxis.set_ticks([], minor=True)
    # Rotate Y Tick Labels for better visibility
    if keep_tick_labels:
        hm.set_yticklabels(hm.get_yticklabels(), rotation=0)
        hm.set_xticklabels(hm.get_xticklabels(), rotation=0)
    else:
        hm.set_yticklabels([])


def _plot_3d_scatter(FRdf, ax, cmap, fig, fill, heatmap_kwargs, index, label):
    """
    Create a scatter plot of excitation (x), inhibition (y), time (z), and firing rate (color).

    Generally a bad visualisation.

    :param FRdf: DataFrame containing firing rate data over time for combinations of E and I synapses.
    :type FRdf: pd.DataFrame
    :param ax: 3D axis to plot upon.
    :type ax: Axes3D
    :param cmap: Colormap to use.
    :type cmap: str
    :param fig: Figure of ax to add colorbar
    :type fig: plt.Figure
    :param fill: Fill the dataframe (True), or keep blanks
    :type fill: bool
    :param heatmap_kwargs: Keywords for sns.heatmap
    :type heatmap_kwargs: dict
    :param index: New index to reindex FRdf
    :type index: pd.Index or pd.MultiIndex
    :param label: Label for firing rate colorbar
    :type label: str
    """
    # resample (dt is 0.025)
    lower_sample = FRdf.iloc[::200, :]
    if fill:
        lower_sample = lower_sample.fillna(method='ffill')
    # change columns to be reference-able
    tt = lower_sample.reindex(columns=index)
    # create a stacked dataframe where 'ex' and 'in' are column names
    stacked = tt.stack(level=['ex', 'in'])
    # 'flatten' the stack (and rename) so that every row is filled for all columns
    df = stacked.reset_index().rename(columns={'level_0': 'time', 0: 'ifr'})
    # color map spec
    vmin, vmax = None, None
    if 'vmax' in heatmap_kwargs.keys():
        vmax = heatmap_kwargs['vmax']
    if 'vmin' in heatmap_kwargs.keys():
        vmin = heatmap_kwargs['vmin']
    colmap = cm.ScalarMappable(norm=colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap)
    colmap.set_array(df[['ifr']])
    # 3D Scatter
    ax.scatter(df[['ex']], df[['in']], df[['time']], marker='.', c=df[['ifr']], s=1)
    fig.colorbar(colmap, ax=ax, label=label)
    ax.set_xlabel('excitation')
    ax.set_ylabel('inhibition')
    ax.set_zlabel('time (ms)')
