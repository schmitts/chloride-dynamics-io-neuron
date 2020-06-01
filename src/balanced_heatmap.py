# coding=utf-8
"""
Run and plot a heatmap for balanced input to find good synapse numbers.
"""
import os
import sys
import logging
from datetime import datetime

import matplotlib.pyplot as plt

from src.config import settings
from src.config.shared import INIT_NEURON
from src.utils.file_io import create_dir, read_ei_file
from src.iocurves.balanced_input.vis import plot_balanced_heatmap
from src.iocurves.balanced_input.synapses_weights import check_num_synapses, check_synapse_strengths

INIT_NEURON()

FOLDER_NAME = "results_balanced_heatmap"

logger = logging.getLogger("balanced input heatmap")


def run_NEURON(exc_weight=1, inh_weight=1, timestamp=False, **kwargs):
    """ Run NEURON simulations where the number of inhibitory and excitatory synapses are varied independently.

    Note that results are added to a txt file, so repeated calls will only run values not present in the file
    already. To force a re-rerun, set timestamp=True

    :param exc_weight: The weighting of the synapses. If `None`, `check_synapse_strengths` is called
    :type exc_weight: int or None
    :param inh_weight: The weighting of the synapses. If `None`, `check_synapse_strengths` is called
    :type inh_weight: int or None
    :param timestamp: Include a timestamp in output file (forces re-run).
    :type timestamp: bool
    :param kwargs: Further keyword arguments to be passed to the sub-method `check_num_synapses`
        - exc_weight, inh_weight
        - exc_syn_list, inh_syn_list
        - exc_freq, inh_freq
        - exc_noise, inh_noise
        - num_runs=5
        - max_num_runs=15
        - upper_bound=15
        - break_bound=True
    :return: List of full-path file names for analysis.
    :rtype: list[str]
    """
    from neuron import h
    date = datetime.strftime(datetime.now(), '%Y-%m-%d_%Hh%M')
    # get output from NEURON (save to file for quicker re-analysis)
    if settings.NEURON_GUI:
        h.showV()
        h.showRunControl()
    create_dir(FOLDER_NAME, timestamp=False)
    try:
        h.hoc_stdout(os.path.join(FOLDER_NAME, f"hoc_stdout_{date}.txt"))
    except RuntimeError:
        pass
    files = []
    for synapse_type in ["distal_KCC2", "proximal_KCC2"]:
        if timestamp:
            f_name = os.path.join(FOLDER_NAME, f"{synapse_type}_{date}.txt")
        else:
            f_name = os.path.join(FOLDER_NAME, f"{synapse_type}.txt")
            if not os.path.exists(f_name):
                with open(f_name, 'w') as _f:
                    _f.write(date)
        logger.info("="*20 + f_name + "="*20)
        with open(f_name, 'r+') as f:
            if exc_weight is None:
                exc_weight, inh_weight = check_synapse_strengths(f)
            logger.info(f"best weights = E:{exc_weight}, I:{inh_weight}")
            check_num_synapses(f, synapse_type, exc_weight, inh_weight, **kwargs)
        files.append(f_name)
    return files


def run_analysis(files=None):
    # Analyse output
    if files is None:
        files = []

    # create a figure per file
    for file in files:
        logger.info(file)
        ei_strength = read_ei_file(file)
        plot_balanced_heatmap(ei_strength, file, detail=True)
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        if sys.argv[1] == "NEURON":
            run_NEURON()
        elif sys.argv[1] == "analysis":
            if len(sys.argv) > 2:
                _files = sys.argv[3:]
            else:
                _files = ["distal_KCC2.txt", "proximal_KCC2.txt"]
            if FOLDER_NAME not in _files[0]:
                _files = [os.path.join(FOLDER_NAME, _f) for _f in _files]
            run_analysis(_files)
        else:
            logger.info("unknown argument")
    else:
        run_analysis(run_NEURON())
