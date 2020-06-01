# coding=utf-8
"""
Run protocols using methods in `iocompare.py`.
The changes in persistence or frequency are specified in `config_stims.hoc`.

This is a holdover from shifting to python from pure hoc scripts.

"""

import numpy as np
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from src.config.shared import INIT_NEURON
from src.utils.nrnconfig import (config_persistent_synapses, single_hz_and_control, config_synapses, balanced_input,
                                 nruns,
                                 )
from src.utils.iocompare import (just_run, compare_diam, compare_pCl, compare_pkcc2, compare_pkcc2_homo,
                                 compare_synapses, compare_weights, compare_cli, compare_pas, compare_dynamic_K,
                                 compare_duration,
                                 )
from src.utils.run_protocol import run_protocol

INIT_NEURON()

change_persist_numbers = None


# This was done in config.hoc and config_stims.hoc instead
#   Instead of None, change_persist_numbers used to be:
#     {"persist_base_config": "persist_base_config*0.1",
#      "ninh_persist_config": "ninh_persist_config*10",
#      "nexc_persist_config": "nexc_persist_config*10"}


# ----------------------------------------------------------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------------------------------------------------------

def protocol_balanced_synaptic_input(distal_balanced_synapses=None, proximal_balanced_synapses=None):
    """
    Runs protocol where synapses have been balanced at 5 Hz (from balanced heatmap).
    For only running balanced frequencies (Inh == Exc), set index_synced: 1. To run multiple Exc Hz per Inh Hz,
    set index_synced to 0.
    """
    config_functions = {balanced_input: {'index_synced': 1}}
    # these were found from balanced_heatmap
    if distal_balanced_synapses is None:
        distal_balanced_synapses = [(200, 30),
                                    (220, 90),
                                    (230, 140),
                                    (240, 200),
                                    (250, 260),
                                    (250, 300),
                                    (260, 400),
                                    (270, 500),
                                    (280, 600),
                                    (290, 700),
                                    (300, 800),
                                    ]
    if proximal_balanced_synapses is None:
        proximal_balanced_synapses = [(260, 30),
                                      (290, 60),
                                      (330, 90),
                                      (370, 120),
                                      (430, 150),
                                      (490, 180),
                                      (530, 210),
                                      (600, 240),
                                      (660, 270),
                                      (720, 300),
                                      (780, 330),
                                      ]
    # not applicable
    skip_if_exists = False

    distal_path = run_protocol(compare_synapses, root="balanced synaptic input distal", timestamp=False,
                               filenames=["distal", "distal_KCC2"],
                               param_list=(skip_if_exists, distal_balanced_synapses),
                               config_functions=config_functions)
    prox_path = run_protocol(compare_synapses, root="balanced synaptic input proximal", timestamp=False,
                             filenames=["proximal", "proximal_KCC2"],
                             param_list=(skip_if_exists, proximal_balanced_synapses),
                             config_functions=config_functions)
    return distal_path, prox_path


def protocol_increasing_synapses_numbers():
    """
    Increases synapses numbers
    """
    config_functions = {
        single_hz_and_control: {"hz": 5, "control": False},
        nruns:                 dict(num_runs=1)
        }
    syn_nums_inh = np.append(0, np.round(np.logspace(1, np.log10(5000), num=10, base=10)).astype(int))
    syn_nums_inh = np.append(syn_nums_inh, np.round(np.geomspace(5000, 10000, num=5)[1:]).astype(int))
    skip = [10, 20, 40, 158]
    for s_val in skip:
        s_idx = np.argmin(np.abs(s_val - syn_nums_inh))
        syn_nums_inh = np.delete(syn_nums_inh, s_idx)
    syn_nums_exc = np.append(0, np.round(np.logspace(1, np.log10(5000), num=19, base=10)).astype(int))
    syn_nums_dict = {
        # 0:    [0, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000],
        # 1:    [0, 10, 20, 50, 100, 200, 500, 1000, 10000],
        # 10:   [0, 10, 20, 50, 100, 200, 500, 1000, 10000],
        # 100:  [0, 10, 20, 50, 100, 200, 500, 1000, 2000, 10000],
        # 1000: [0, 100, 200, 500, 1000, 2000, 5000, 10000],
        # 2000: [0, 100, 200, 500, 1000, 2000, 5000, 10000]
        }
    for i in np.append(syn_nums_inh, np.round(np.geomspace(5000, 10000, num=5)[1:]).astype(int)):
        syn_nums_dict[i] = np.append(syn_nums_exc, np.round(np.geomspace(1000, 100000, num=15)).astype(int))
    syn_nums = []
    for i, e_list in syn_nums_dict.items():
        for e in e_list:
            syn_nums.append((e, i))
    # syn_nums = [s for s in itertools.product(syn_nums_exc, syn_nums_inh)]

    run_protocol(compare_synapses, root="Increasing synapse numbers", timestamp=False,
                 filenames=["distal", "distal_KCC2"],
                 param_list=(True, syn_nums), config_functions=config_functions)

    for i in np.round(np.geomspace(630, 2507, num=7)).astype(int):
        syn_nums_dict[i] = np.append(0, np.round(np.geomspace(1000, 100000, num=15)).astype(int))
    syn_nums = []
    for i, e_list in syn_nums_dict.items():
        for e in e_list:
            syn_nums.append((e, i))
    # run_protocol(compare_synapses, root="Hz", name="Increasing synapse numbers", timestamp=False,
    #          filenames=["proximal", "proximal_KCC2"],
    #          param_list=(True, syn_nums), config_functions=config_functions)


def protocol_increasing_weights():
    """
    Increases weights of synapse
    """
    config_functions = {
        single_hz_and_control: {"hz": 5, "control": False},
        config_synapses:       {"exc": 100, "inh": 100}
        }
    exc_weights = np.array([0.1, 1, 10])
    exc_base_weight = 1
    inh_weights = np.array([0.1, 1, 10])
    inh_base_weight = 1
    # set the actual weights to be as a multiple of the base weights
    weights_dict = {
        'exc_weights': (exc_weights*exc_base_weight).tolist(),
        'inh_weights': (inh_weights*inh_base_weight).tolist()
        }
    run_protocol(compare_weights, root=protocol_increasing_weights.__name__, timestamp=False,
                 filenames=["proximal", "proximal_KCC2", "distal", "distal_KCC2"], param_list=[weights_dict],
                 config_functions=config_functions)


def protocol_increasing_persistence():
    """
    Increase the conductance of persistent synapses.
    `just_run` doesn't tak any parameters.
    Values are configured in `config_stims.hoc`
    """
    config_functions = {config_persistent_synapses: change_persist_numbers}
    run_protocol(just_run, root="persistent", timestamp=False, filenames=["proximal", "proximal_KCC2"],
                 param_list=None,
                 config_functions=config_functions)
    run_protocol(just_run, root="persistent", timestamp=True, filenames=["distal", "distal_KCC2"],
                 param_list=None,
                 config_functions=config_functions)


def protocol_increasing_persistence_change_diam():
    """
    Change diameter of persistent synapse model

    """
    config_functions = {config_persistent_synapses: change_persist_numbers}
    run_protocol(compare_diam, root="persistent_diam", timestamp=False, filenames=["distal", "distal_KCC2"],
                 param_list=[0.5, 1.0, 1.5, 2.0], config_functions=config_functions)


def protocol_increasing_persistence_change_pcl():
    """
    Change base chloride of persistent synapse model

    """
    config_functions = {config_persistent_synapses: change_persist_numbers}
    run_protocol(compare_pCl, root="persistent_pcl", timestamp=False, filenames=["distal", "distal_KCC2"],
                 param_list=None,
                 config_functions=config_functions)


def protocol_increasing_persistence_change_pkcc2():
    """
    Change KCC2 strength of persistent synapse model

    """
    base = 1.9297e-5
    pkcc2_list = [base/2, base, base*2, base*4]
    config_functions = {config_persistent_synapses: change_persist_numbers}
    run_protocol(compare_pkcc2, root="persistent_pkcc2", timestamp=False, filenames=["distal_KCC2"],
                 param_list=[pkcc2_list[-1]], config_functions=config_functions)


def protocol_increasing_persistence_change_pkcc2_homo():
    """
    Change homogeneity of KCC2 strength of persistent synapse model
    """
    config_functions = {config_persistent_synapses: change_persist_numbers}
    run_protocol(compare_pkcc2_homo, root="persistent_pkcc2_homo", timestamp=False, filenames=["distal", "distal_KCC2"],
                 param_list=None, config_functions=config_functions)


def protocol_increasing_persistence_change_pas():
    """
    Change leak conductance (input resistance) of persistent synapse model
    """
    config_functions = {config_persistent_synapses: change_persist_numbers}
    run_protocol(compare_pas, root="persistent_pas", timestamp=False, filenames=["distal", "distal_KCC2"],
                 param_list=None,
                 config_functions=config_functions)


def protocol_increasing_persistence_change_duration():
    """
    Change duration of persistent synapse model
    """
    config_functions = {config_persistent_synapses: change_persist_numbers}
    run_protocol(compare_duration, root="persistent_duration", timestamp=False, filenames=["distal"],
                 param_list=[10000.0],
                 config_functions=config_functions)


def protocol_increasing_persistence_change_dynamic_K():
    """
    Change duration of persistent synapse model
    """
    config_functions = {config_persistent_synapses: change_persist_numbers}
    run_protocol(compare_dynamic_K, root="persistent_K", timestamp=False, filenames=["distal"],
                 param_list=[True, False],
                 config_functions=config_functions)


def protocol_increasing_persistence_change_cli():
    """
    Change internal chloride of persistent synapse model
    (only works for non-KCC2 files)

    """
    config_functions = {config_persistent_synapses: change_persist_numbers}
    run_protocol(compare_cli, root="persistent_cli", timestamp=False, filenames=["distal", "proximal"], param_list=None,
                 config_functions=config_functions)


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        if sys.argv[1] == "balanced":
            protocol_balanced_synaptic_input()
        elif sys.argv[1] == "persistence":
            if len(sys.argv) > 2:
                choose = {
                    "pkcc2_homo": protocol_increasing_persistence_change_pkcc2_homo,
                    "duration":   protocol_increasing_persistence_change_duration,
                    "dynamic_K":  protocol_increasing_persistence_change_dynamic_K,
                    "pas":        protocol_increasing_persistence_change_pas,
                    "diam":       protocol_increasing_persistence_change_diam,
                    "pkcc2":      protocol_increasing_persistence_change_pkcc2,
                    "pcl":        protocol_increasing_persistence_change_pcl,
                    "cli":        protocol_increasing_persistence_change_cli,
                    }
                for choice in sys.argv[2:]:
                    print("choice = {}".format(choice))
                    choose[choice]()
            else:
                print("protocol_increasing_persistence")
                protocol_increasing_persistence()
        elif sys.argv[1] == "synapses":
            protocol_increasing_synapses_numbers()
        else:
            print("please choose one of 'balanced' 'persistence' 'synapses'")
    else:
        print("please choose one of 'balanced' 'persistence' 'synapses'")
