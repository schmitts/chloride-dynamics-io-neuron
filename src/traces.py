# coding=utf-8
"""
Generate and plot voltage and/or chloride ion concentration trace for neuron(s) using hoc_files.
"""
import logging
import pandas as pd
import seaborn as sns
import numpy as np
import scipy.interpolate as interpolate
from matplotlib import pyplot as plt
from neuron import h
from numbers import Number

from src.config import settings, shared
from src.config.settings import COLOR
from src.utils.file_io import create_dir
from src.utils.nrnpy import get_base_vm_cli, load_file
from src.utils.timing import current_time

logger = logging.getLogger("trace")

shared.INIT_NEURON()


def get_trace(file_name, synapse_type=0, synapse_numbers=(10, 10),
              hz=None, space_plot=False, vm=-65.0, cli=5.0, ldend_diam=0.5):
    """
    Get the voltage trace and cli trace from running a hoc file.

    Optionally, specify space_plot to get the cli along the entire neuron (nseg dependent).

    `vm` and `cli` should be specified by running `get_base_vm_cli` beforehand and passing in those values.

    :param file_name: hoc-specified neuron to load (from hoc_files/cells)
    :type file_name: str
    :param synapse_type: 0 if 'f-in' NETSTIM synapses. 1 if 'persistetnt' fluctuating conductance synapses.
    :type synapse_type: int
    :param synapse_numbers: Number of (E,I) synapses.
    :type synapse_numbers: tuple of int
    :param hz: Input to the synapses. Keys used are 'in' and 'ex'.
        Either frequency-based (`synapse_type` is `0`) or relative conductance (`synapse_type` is `1`).
    :type hz: dict of str
    :param space_plot: Compute cli at every location.
    :type space_plot: bool
    :param vm: Initial membrane potential
    :type vm: float
    :param cli: Initial chloride ion concentration
    :type cli: float
    :param ldend_diam: Diameter of distal dendrite.
    :type ldend_diam: float
    :return: Time array, voltage array, chloride array,
        space plot dict data of form `{distance from soma: chloride array}`
    :rtype: tuple[h.Vector, h.Vector, h.Vector, dict of float:h.Vector]
    """
    if hz is None:
        hz = {'in': 5, 'ex': 5}
    load_file(file_name)

    h.changeSynapseType(synapse_type)
    h.newSynapses(synapse_numbers[0], synapse_numbers[1])
    if synapse_type == 0:
        h.inPy(hz['in'], 1, 1)
        h.ex(hz['ex'], 1, 1)
    else:
        h.inPy(hz['in'])
        h.ex(hz['ex'])
    if "distal" in file_name:
        compartment = h.ldend
    elif "proximal" in file_name:
        compartment = h.bdend
    else:
        compartment = h.soma

    h("forall {" +
      "cli={} ".format(cli) +
      "cli0_KCC2={}".format(cli) +
      "}")
    for sec in h.allsec:
        for seg in sec:
            seg.cli = cli
    h.cli0_cl_ion = cli
    h.ldend.diam = ldend_diam
    # sort on key
    data = {}
    if space_plot:
        h.distance(0, 1, sec=h.soma)
        for sec in h.allsec:
            for seg in sec:
                # label = f"{sec.name()}({seg.x:.5f})"
                label = h.distance(seg.x, sec=sec)
                if 'dend' in sec.name():
                    label *= -1
                data[label] = h.Vector()
                data[label].record(sec(seg.x)._ref_cli)

    h.tstop = 1000.0
    time_past = current_time('ms')
    logger.info("running {}...".format(file_name))
    h.v_init = vm
    h.finitialize(vm)
    t_vec = h.Vector()
    t_vec.record(h._ref_t)
    v_vec = h.Vector()
    v_vec.record(h.axon(0.5)._ref_v)
    cli_vec = h.Vector()
    cli_vec.record(compartment(0.5)._ref_cli)

    h.run()

    logger.info("time taken: {}ms".format(current_time('ms') - time_past))
    logger.info("no. spikes = {}".format(h.apc.n))

    return t_vec, v_vec, cli_vec, data


def space_from_dict(data: dict):
    """Convert `data` dictionary of {x:y} to x, y arrays sorted by x."""
    x, y = [], []
    for key, val in sorted(data.items()):
        x.append(key)
        y.append(val.mean())
    return np.array(x), np.array(y)


def smooth(x, y, N):
    """Smoothing function by interpolating x and y to N points.

    :param x: x-coordinates
    :type x: list or np.ndarray
    :param y: y-coordinates (same size as x)
    :type y: list or np.ndarray
    :param N: Size of new x, y arrays
    :type N: int
    :return: new x and y arrays
    :rtype: tuple[np.ndarray]
    """
    x_sm = np.array(x)
    y_sm = np.array(y)

    x_smooth = np.linspace(x_sm.min(), x_sm.max(), N)
    t, c, k = interpolate.splrep(x, y, s=0, k=4)
    spline = interpolate.BSpline(t, c, k, extrapolate=False)
    return x_smooth, spline(x_smooth)


def space_plot(data, ax=None, **kwargs):
    """Create 2D space plot of a neuron with a singular color."""
    x, y = space_from_dict(data)
    x_sm, y_sm = smooth(x, y, 200)
    ax.plot(x_sm, y_sm, **kwargs)
    ax.set(ylabel=f"{settings.CLI} ($mM$)", xlabel='Distance from soma ($\mu m$)')


def plot_v_trace(*args, title="", same_ax=False, show_cli=False, cmap=None, mean=False, dt=20.,
                 time_points=(20., 100., 500., 1000.),
                 **kwargs):
    """
    Plot voltage from passed data (in `*args`) with optional cli plot.

    :param args: List of data dictionaries.
    :type args: list[dict]
    :param title: Figure title
    :type title: str
    :param same_ax: Plot on the same axis (default: False)
    :type same_ax: bool
    :param show_cli:
    :type show_cli: bool
    :param cmap: Colormap to use.
    :type cmap: sns.palettes._ColorPalette or list
    :param mean: Plot the average of the voltage and [Cl-]i too (default: False).
    :type mean: bool
    :param dt: Time step to underline for each time point in time_points
    :type dt: float
    :param time_points:
    :type time_points: Iterable[float]
    :param kwargs: further arguments to pass to `Axes.plot`
    :type kwargs: dict
    :return: Figure, Axes array pair made.
    :rtype: (matplotlib.figure.Figure, np.ndarray of matplotlib.axes.Axes)
    """
    # Voltage plot

    if mean: same_ax = True
    nrows = 1 if same_ax else len(args)
    fig, ax = plt.subplots(nrows=nrows + show_cli, ncols=1, figsize=(settings.PAGE_W_half, settings.PAGE_H_3rd),
                           sharex=True, squeeze=False, gridspec_kw={"height_ratios": [1]*nrows + [nrows]*show_cli})
    fig.suptitle(title)
    y = -75
    cli_ax = None
    vm_avg = []
    cli_avg = []
    for i, data in enumerate(args):
        t_vec, v_vec, cli_vec, space_data = data
        vm_avg.append(np.array(v_vec))
        axis = ax[0, 0] if same_ax else ax[i, 0]
        c = None if cmap is None else cmap[i]
        if c is None and mean:
            c = 'k'
        alpha = 0.2 if mean else 1
        axis.plot(t_vec, v_vec, c=c, lw=0.75, alpha=alpha, **kwargs)
        axis.set_xlim([0, h.tstop])
        # plot specific points on the voltage plot (for heatmap trace)
        if not same_ax or (same_ax and i == 0):
            for tp in time_points:
                axis.axhline(y, xmin=(tp - dt)/1000.0, xmax=tp/1000.0, color='k', lw=0.5)
        if show_cli:
            cli_ax = ax[-1, 0]
            cli_ax.plot(t_vec, cli_vec, c=c, alpha=alpha, **kwargs)
            cli_avg.append(np.array(cli_vec))
    if mean:
        mean_v = np.mean(vm_avg, axis=0)
        axis.plot(t_vec, mean_v, c='k', lw=0.75, alpha=1, **kwargs)
        if show_cli:
            mean_cli = np.mean(cli_avg, axis=0)
            cli_ax.plot(t_vec, mean_cli, c='k', lw=0.75, alpha=1, **kwargs)
    ax[-1, 0].set_xlabel("Time (ms)")
    if show_cli:
        ax[(len(ax) - 1)//2, 0].set_ylabel("Membrane\npotential\n(mV)")
        ax[-1, 0].set_ylabel("Chloride\nconcentration\n(mM)")
    else:
        ax[-1, 0].set_ylabel("Membrane\npotential\n(mV)")
    return fig, ax


def figure_cli_distribution(syn_n_dist=(250, 300), syn_n_prox=(330, 90), savefig=False):
    """
    Plot the distribution of [Cl-]i along a neuron for distal and proximal inhibitory input.

    :param syn_n_dist: Number of (E, I) synapses for the distal input simulation.
    :type syn_n_dist: tuple
    :param syn_n_prox: Number of (E, I) synapses for the proximal input simulation.
    :type syn_n_prox: tuple
    :param savefig: Whether to save the figure to results_plots.
    :type savefig: bool
    """
    txt = f"figure_cli_distribution(syn_n_dist={syn_n_dist}, syn_n_prox={syn_n_prox})"
    logger.info(txt)
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 3.5), sharey=True, sharex=True)
    fig.subplots_adjust(hspace=0.5)

    for filename, syn_n, _ax in zip(['distal_KCC2', 'proximal_KCC2'], [syn_n_dist, syn_n_prox], ax):
        # steady state vm and [Cl-]i
        vm, cli = get_base_vm_cli(filename, load=True)
        trace_kwargs = dict(synapse_type=0,  # frequency-based
                            synapse_numbers=syn_n,  # (E,I)
                            space_plot=True)  # compute spatial component for cli

        # run simulations for different frequencies
        _, _, _, inh = get_trace(filename, vm=vm, cli=cli, hz={'in': 5, 'ex': 0}, **trace_kwargs)
        _, _, _, exc = get_trace(filename, vm=vm, cli=cli, hz={'in': 0, 'ex': 5}, **trace_kwargs)
        _, _, _, equal = get_trace(filename, vm=vm, cli=cli, hz={'in': 5, 'ex': 5}, **trace_kwargs)

        logger.info(f"space plot {filename}")
        space_plot(exc, ax=_ax, color=COLOR.E, label='5 : 0')
        space_plot(inh, ax=_ax, color=COLOR.I, label='0 : 5')
        space_plot(equal, ax=_ax, color=COLOR.E_I, label='5 : 5')
        # display some info on figure
        name = filename.split("_")[0].capitalize()
        _ax.set_title(f"{name} input (E: I) (# of synapses)\n"
                      f"{syn_n[0]}:{syn_n[1]:>3.0f}",
                      va='top', fontsize='medium')
    # light adjustments
    ax[0].set_xlabel("")
    # ax[0].set_xticklabels([])
    ax[0].set_xlim(-550, 500)  # 0 is soma(0)
    ax[-1].legend(title='Input (E : I) (Hz)', loc='upper right', frameon=False)
    if savefig:
        create_dir("results_plots", timestamp=False)
        fig.savefig(f"results_plots/{txt}.png", bbox_inches='tight', frameon=False)
        fig.savefig(f"results_plots/{txt}.pdf", bbox_inches='tight', frameon=False)


def figure_cli_heatmaps(distal, proximal, freqs=(5, 10, 25, 50),
                        n_trials=1, show_fr=False, vmin=0, vmax=None, savefig=False):
    """
    Heatmaps of [Cl-]i for different frequencies and (E,I) synapse pairs and for both proximal and distal.

    :param distal: List (E,I) synapse number pairs for distal input.
    :type distal: list of (tuple of (float))
    :param proximal: List (E,I) synapse number pairs for proximal input.
    :type proximal: list of (tuple of (float))
    :param freqs:
    :type freqs: list of Number or tuple of Number
    :param n_trials: Nomber of repeated runs.
    :type n_trials: int
    :param show_fr: Show the firing rate to the right of the heatmaps when plotting.
        True shows the firing rate as an arrow and text.
        Any value greater than 1 shows the firing rate as a heatmap cell (different colormap)
    :type show_fr: bool or int
    :param vmin: Minimum [Cl-]i for color range.
    :type vmin: float
    :param vmax: Maximum [Cl-]i for color range.
    :type vmax: float
    :param savefig: Whether to save the figure to results_plots
    :type savefig: bool
    """
    logger.info(f"figure_cli_heatmaps(distal={distal}, proximal={proximal}, "
                f"n_trials={n_trials}, show_fr={show_fr})")

    fig, ax2d = plt.subplots(nrows=2, ncols=len(distal) + 1 + int(show_fr > 1), figsize=(8, 5),
                             gridspec_kw={'width_ratios': [15]*len(distal) + [2]*(1 + int(show_fr > 1))})
    cmap = sns.cubehelix_palette(16, start=2.7, rot=-.2, light=0.98, dark=0.40, as_cmap=True)
    cmap_fr = sns.color_palette("Reds" if show_fr > 1 else "magma", n_colors=200, desat=1)
    vmax_fr = 5
    # h.hoc_stdout("hoc_output_traces.txt")
    for fdx, (filename, synapse_numbers) in enumerate(zip(["distal_KCC2", "proximal_KCC2"],
                                                          [distal, proximal])):
        ax = ax2d[fdx, :]
        cbar_ax = ax[-1 - int(show_fr > 1)]
        vm, cli = get_base_vm_cli(filename, load=True)

        d_cli = {}
        fr = {}
        stddev = {}
        for i, syn_n in enumerate(synapse_numbers):
            kwargs = dict(synapse_type=0,
                          synapse_numbers=syn_n,
                          space_plot=True)
            for hz in freqs:
                n_spikes = []
                sec_means = [0, 0, 0, 0]  # 4 sections [ldend, bdend, soma, axon]
                for t in range(n_trials):
                    logger.info(f"filename={filename} syn_n={syn_n} hz={hz} t={t}")
                    _, _, _, data = get_trace(filename, vm=vm, cli=cli, hz={'in': hz, 'ex': hz}, **kwargs)
                    n_spikes.append(h.apc.n)
                    x, y = space_from_dict(data)
                    # separate x into regions based on known lengths
                    ldend = x <= -50
                    bdend = np.logical_and(-50 <= x, x <= 0)
                    soma = np.logical_and(0 <= x, x <= 15)
                    axon = x >= 15
                    for s, sec in enumerate([ldend, bdend, soma, axon]):
                        sec_means[s] += y[sec].mean()
                sec_means = [m/n_trials for m in sec_means]
                d_cli[(f"{syn_n[0]:>3.0f}:{syn_n[1]:>3.0f}", hz)] = sec_means
                fr[(f"{syn_n[0]:>3.0f}:{syn_n[1]:>3.0f}", hz)] = sum(n_spikes)/n_trials
                stddev[(f"{syn_n[0]:>3.0f}:{syn_n[1]:>3.0f}", hz)] = np.std(n_spikes)

        df = pd.DataFrame.from_dict(d_cli, orient='index',
                                    columns=["Distal\nDendrite", "Proximal\nDendrite", "Soma", "Axon"])
        df_fr = pd.DataFrame.from_dict(fr, orient='index',
                                       columns=["Output"])
        df = df.reindex(pd.MultiIndex.from_tuples(df.index))
        df_fr = df_fr.reindex(pd.MultiIndex.from_tuples(df_fr.index))
        vmin = vmin or df.values.min()
        vmax = vmax or df.values.max()
        vmin_fr = 0
        vmax_fr = max(vmax_fr, df_fr.values.max())*1.1  # give a 10% buffer for the cmap

        logger.info("plotting cli_heatmaps")
        for i, syn_n in enumerate(df.index.levels[0]):
            df_syn = df.loc[syn_n]
            df_fr_syn = df_fr.loc[syn_n]
            sns.heatmap(df_syn, ax=ax[i], annot=False, fmt=".1f", square=True,
                        annot_kws=dict(fontsize='xx-small'),
                        vmin=vmin, vmax=vmax,
                        cmap=cmap,
                        cbar=(i == 0),
                        cbar_ax=None if i > 0 else cbar_ax,
                        )

            if show_fr > 1:
                new_df = pd.concat([df_syn, df_fr_syn], axis='columns')
                mask = np.ones(df_syn.shape)
                mask = np.append(mask, np.zeros(df_fr_syn.shape), axis=1)
                cbar_fr_ax = ax[-1]
                sns.heatmap(new_df, ax=ax[i], annot=True, fmt=".1f", square=False,
                            annot_kws=dict(fontsize='xx-small'),
                            vmin=vmin_fr, vmax=vmax_fr,
                            cmap=cmap_fr,
                            mask=mask,
                            cbar=(i == 0),
                            cbar_ax=cbar_fr_ax,
                            cbar_kws={'label': "Firing rate (Hz)"},
                            )
            elif show_fr:
                for j, _fr in enumerate(df_fr_syn['Output']):
                    idx = (len(cmap_fr) - 1)*(_fr - vmin_fr)//(vmax_fr - vmin_fr)
                    c = cmap_fr[int(idx)]
                    text = ax[i].annotate(f"{_fr:>2.1f}", xy=(4, j + 0.5), xytext=(4.6, j + 0.5), color=c, alpha=1,
                                          arrowprops={'arrowstyle': '<-'}, fontsize='x-small', va='center')
                    # path effects can sometimes make text clearer, but it can also make things worse...
                    # text.set_path_effects([path_effects.Stroke(linewidth=0.5, foreground='black', alpha=0.5),
                    #                        path_effects.Normal()])
            ax[i].set_xticklabels(ax[i].get_xticklabels(), fontsize="small", rotation=45)
            ax[i].set_title(syn_n, fontsize='small')
            if i == 0:
                ax[i].set(ylabel='Balanced Input (Hz)')
                ax[i].set_yticklabels(df_syn.index, rotation=0)
            else:
                ax[i].set_yticklabels([])
            if i == 0:
                cbar_ax.set_xlabel(f'{settings.CLI} (mM)', ha='center')
                cbar_ax.xaxis.set_label_position('top')
    fig.tight_layout()
    if savefig:
        create_dir("results_plots", timestamp=False)
        fig.savefig(f"results_plots/figure_cli_heatmaps_{int(show_fr)}.png", bbox_inches='tight', facecolor='None')
        fig.savefig(f"results_plots/figure_cli_heatmaps_{int(show_fr)}.pdf", bbox_inches='tight', facecolor='None')


def figure_v_traces(inh_region="distal", KCC2=True, show_cli=False, base=True,
                    synapse_type=0, synapse_numbers=(100, 100), hz=None, mean=False,
                    savefig=False, **kwargs):
    """
    Display voltage traces for an inhibitory input region (e.g. "distal") and (E, I) synapse pairs.

    Optionally include an axis for [Cl-]i.

    :param inh_region: Region for inhibitory synapses. Only those in hoc_files/cells/ are supported.
    :type inh_region: str
    :param KCC2: Open the '_KCC2' version of the hoc file instead if True.
    :type KCC2: bool
    :param show_cli: Include a plot for [Cl-]i. Only applicable when KCC2 is True.
    :type show_cli: bool
    :param base: Find steady-state values for vm and cli (`True`),
        provide the steady-state values (`(<vm value>, <cli value>)`),
        use defaults of -71 mV and 4.25 mM for vm and cli, respectively(`False`).
    :type base: bool or tuple of float
    :param synapse_type: Use frequency-based synapses (`0`) or persistant synapses (`1`).
    :type synapse_type: int
    :param synapse_numbers: Synapses numbers for the neuron in (E, I) format. A list of (E,I) pairs can be provided
        for multiple traces.
    :type synapse_numbers: list of (tuple of (int)) or tuple of int
    :param hz: Frequency of synapses (if `synapse_type` is `0`) or the relative conductance (if `synapse_type` is `1`).
    :type hz: dict
    :param mean: Additionally plot the mean vm and cli values using `plot_v_trace`
    :type mean: bool
    :param savefig: Save the figure to results_plots
    :type savefig: bool
    :param kwargs: Keyword arguments to pass to `get_trace`
    :type kwargs: dict
    :return: Steady-state voltage and [Cl-]i used for this simulation
    :rtype:
    """
    if hz is None:
        hz = {'in': 5, 'ex': 5}
    logger.info(f"figure_v_traces(inh_region={inh_region} KCC2={KCC2}, synapse_numbers={synapse_numbers}, hz={hz})")
    cmap_name = "Blues" if inh_region == "distal" else "Greens"
    cmap = sns.color_palette(cmap_name, n_colors=len(synapse_numbers))
    if mean: cmap = None
    filename = inh_region
    if KCC2:
        filename += "_KCC2"
    else:
        show_cli = False
    logger.info("getting base vm cli")
    if base:
        if type(base) is tuple:
            vm, cli = base
        else:
            vm, cli = get_base_vm_cli(f"{inh_region}_KCC2", load=True)
    else:
        vm, cli = -71., 4.25

    if type(synapse_numbers) is tuple:
        synapse_numbers = [synapse_numbers]

    dynamic_data = []
    spikes = []
    for syn_num in synapse_numbers:
        trace_kwargs = dict(synapse_type=synapse_type,
                            synapse_numbers=syn_num,
                            hz=hz,
                            space_plot=False, **kwargs)
        dynamic_data.append(get_trace(filename, vm=vm, cli=cli, **trace_kwargs))
        spikes.append(h.apc.n)
    logger.info("plotting voltage")
    fig, ax = plot_v_trace(*dynamic_data, show_cli=show_cli, cmap=cmap, mean=mean)
    if not mean:
        for i, (spike_num, syn_num) in enumerate(zip(spikes, synapse_numbers)):
            ax[i, 0].text(ax[i, 0].get_xlim()[1], 0, f"{spike_num:.0f}", ha="left")
            ax[i, 0].set_title(str(syn_num).replace("(", "").replace(")", ""),
                               fontsize='small', va='top')
    title = filename + str(hz)
    title = title.replace("{", "\n(").replace("}", ")").replace("'", "").replace(": ", "=")
    if synapse_type == 0:
        title = title.replace(",", " Hz,").replace(")", " Hz)")
    fig.suptitle(title, fontsize='medium')
    fig.align_ylabels(ax[:, 0])
    sns.despine(fig, bottom=True, left=True)
    for _ax in ax[:-1, 0]:
        _ax.xaxis.set_visible(False)
        _ax.yaxis.set_visible(False)
    sns.despine(ax=ax[-1, 0])
    fig.subplots_adjust(left=0.15, bottom=0.15)
    if savefig:
        create_dir("results_plots")
        fig.savefig(f"results_plots/figure_trace_{title}.png", bbox_inches='tight', facecolor='None')
        fig.savefig(f"results_plots/figure_trace_{title}.pdf", bbox_inches='tight', facecolor='None')
    return vm, cli


if __name__ == "__main__":
    # options
    _distal_balanced_synapses = [(200, 30), (220, 90), (230, 140), (240, 200), (250, 260), (250, 300), (260, 400),
                                 (270, 500), (280, 600), (290, 700), (300, 800)]
    _proximal_balanced_synapses = [(260, 30), (290, 60), (330, 90), (370, 120), (430, 150), (490, 180), (530, 210),
                                   (600, 240), (660, 270), (720, 300), (780, 330)]
    # chosen
    distal_balanced_synapses = [(200, 30), (230, 140), (250, 300), (270, 500), (300, 800)]
    proximal_balanced_synapses = [(260, 30), (290, 60), (330, 90), (490, 180), (780, 330)]
    mid_d = len(distal_balanced_synapses)//2
    mid_p = len(proximal_balanced_synapses)//2
    distal = [distal_balanced_synapses[0], distal_balanced_synapses[mid_d], distal_balanced_synapses[-1]]
    proximal = [proximal_balanced_synapses[0], proximal_balanced_synapses[mid_p], proximal_balanced_synapses[-1]]

    # figure 2 shape plot for cli along neuron
    figure_cli_distribution(distal[1], proximal[1])
    # figure 2 heatmaps for cli concentration at different input frequencies
    figure_cli_heatmaps(distal, proximal, show_fr=True)

    # trace for fig 5 heatmap
    vm, cli = figure_v_traces(inh_region="distal", KCC2=True, show_cli=True, synapse_type=1,
                              synapse_numbers=[(100, 100)]*10, hz={'in': 4, 'ex': 12},
                              mean=False)

    figure_v_traces(inh_region="distal", synapse_numbers=distal,
                    hz={'in': 5, 'ex': 5}, show_cli=True)
    figure_v_traces(inh_region="proximal", synapse_numbers=proximal,
                    hz={'in': 5, 'ex': 5}, show_cli=True)
    figure_v_traces(inh_region="distal", synapse_numbers=distal_balanced_synapses,
                    hz={'in': 20, 'ex': 20}, show_cli=True)
    figure_v_traces(inh_region="proximal", synapse_numbers=proximal_balanced_synapses,
                    hz={'in': 20, 'ex': 20}, show_cli=True)
    figure_v_traces(inh_region="distal", KCC2=False, synapse_numbers=distal_balanced_synapses,
                    hz={'in': 20, 'ex': 20})
    figure_v_traces(inh_region="proximal", KCC2=False, synapse_numbers=proximal_balanced_synapses,
                    hz={'in': 20, 'ex': 20})

    logger.info("showing plots")
    plt.show()
