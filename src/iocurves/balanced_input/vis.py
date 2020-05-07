"""
Visualise data from balanced input-ouput
"""
import logging
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

logger = logging.getLogger("balanced input vis")


def cond_change_to_max(spike, max_spike):
    if spike > max_spike:
        return 10
    else:
        return spike


def plot_balanced_heatmap(constants, file_name, detail=False):
    import pandas as pd
    import seaborn as sns
    cmap_bwr = sns.diverging_palette(220, 20, sep=10, n=100, as_cmap=True)
    for (key, spikes) in constants.items():
        logger.info(f"weights = {key}")
        spike_list, stddev_list, n_list = zip(*spikes.values())
        n_arr = np.array(n_list)
        if detail:
            # heatmapkw = dict(annot=True, fmt="1.2f", annot_kws={"size": 7})
            heatmapkw = dict()
            file_name = file_name.replace(".txt", "_detail.txt")
            figsize = (8, 3.5)
        else:
            heatmapkw = dict()
            figsize = (4, 3.5)
        spike_list = [cond_change_to_max(spike, 10) for spike in spike_list]

        fig, ax = plt.subplots(1, 1 + int(detail), figsize=figsize,
                               sharex=True, sharey=True, squeeze=False)

        for i, (label, lst, cmap, ax_plot) in enumerate(zip(["Firing Rate (Hz)", "std dev"],
                                                            [spike_list, stddev_list],
                                                            [cmap_bwr, 'plasma'],
                                                            ax[0, :])):
            ser = pd.Series(lst, index=pd.MultiIndex.from_tuples(spikes.keys()))
            # if we haven't really sampled the data (high firing rate), set to nan value
            ser.iloc[n_arr == 0] = np.nan
            df = ser.unstack()
            # color missing data
            sns.heatmap(df.isnull(), cmap=sns.light_palette('lavender', 1), square=True, ax=ax_plot, cbar=False,
                        xticklabels='auto', yticklabels='auto', lw=0.01, **heatmapkw)
            # color actual data
            hm = sns.heatmap(df, cmap=cmap, square=True, ax=ax_plot, cbar_kws={'label': label},
                             xticklabels='auto', yticklabels='auto', mask=df.isnull(), **heatmapkw)

            ax_plot.invert_yaxis()
            # ax_plot.set_yticklabels(ax_plot.get_yticklabels(), rotation=45)
            ax_plot.set_xticklabels(ax_plot.get_xticklabels(), rotation=60)
            ax_plot.set_xlabel("Number of $GABA_A$ synapses")
            ax_plot.set_ylabel("Number of mixed $NMDA/AMPA$ synapses")

        fig.suptitle("file:" + file_name, fontsize='small')

        fig.tight_layout()
        for format in ["eps", "png"]:
            save_name = file_name.replace(".txt", f".{format}")
            fig.savefig(save_name)

        keys = list(spikes.keys())
        closest = (abs(5 - spikes[keys[0]][0]), keys[0][0], keys[0][1])
        inh_closest = {}
        logger.info("exc \t inh \t spike \t stddev \t number of trials")
        for ((exc, inh), (spike, stddev, num)) in spikes.items():
            _dif = abs(5 - spike)
            if _dif < closest[0]:
                closest = (_dif, exc, inh)
            if _dif < 1:
                logger.info(f"{exc:>3.0f} \t {inh:>3.0f} \t {spike:>5.2f} \t {stddev:>4.2f} \t {num:2.0f})")
            if inh not in inh_closest or _dif < inh_closest[inh][0]:
                inh_closest[inh] = (_dif, exc)

        logger.info("\t".join(list(zip(["offset:", "excitation", "inhibition"], closest))))

        for inh, (_dif, exc) in inh_closest.items():
            logger.info("({},{}) = {}".format(exc, inh, _dif))


def plot_balanced_iocurves(df: pd.DataFrame, cmap='Blues', synapse_list=None):
    sns.set_context('talk')
    logger.info(f"plotting {synapse_list if synapse_list is not None else ''}")
    # create balanced column
    df['Balanced Input (Hz)'] = df[['E (Hz)', 'I (Hz)']].apply(lambda x: x[0] if x[0] == x[1] else np.nan, axis=1)
    # and remove unbalanced entries
    df = df[~df['Balanced Input (Hz)'].isnull()]

    # create joint synapse column
    df['Synapses (E:I)'] = df[['E (Synapses)', 'I (Synapses)']]\
        .apply(lambda x: ':'.join([f"{_x:>3.0f}" for _x in x]), axis='columns')
    # and filter by passed `synapse_list`
    if synapse_list is not None:
        if type(synapse_list[0]) is tuple:
            # convert from numbers to str for finding in 'Synapses (E:I)' and hue order
            synapse_list = [f"{_e:>3.0f}:{_i:>3.0f}" for _e, _i in synapse_list]
        df = df[df['Synapses (E:I)'].isin(synapse_list)]

    # create a line plot with stddev error bands
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.lineplot(x='Balanced Input (Hz)', y='Firing Rate (Hz)',
                 hue='Synapses (E:I)', hue_order=synapse_list,
                 style='KCC2', style_order=[True, False],
                 palette=cmap, ci='sd',
                 data=df)
    leg = ax.legend(frameon=True, facecolor='w', edgecolor='None')

    xeqy = np.arange(ax.get_xlim()[1])
    ax.plot(xeqy, xeqy, color='k', alpha=0.5, lw=0.5, zorder=99)
    ax.set_xlim(0, 60)
    ax.set_ylim(0)
    sns.despine(ax=ax)
    return ax
