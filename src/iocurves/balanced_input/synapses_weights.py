"""
Simulation methods for running balanced input (Hz) for different excitation+inhibition synapses numbers, and strength
 of synapses (aka weights)
"""
import logging
import numpy as np
from neuron import h

from src.utils.nrnpy import load_file

logger = logging.getLogger("balanced input")


def check_num_synapses(f, synapse_type, exc_weight, inh_weight,
                       exc_syn_list=None, inh_syn_list=None,
                       exc_freq=5, inh_freq=5,
                       exc_noise=1, inh_noise=1,
                       num_runs=5, max_num_runs=15,
                       upper_bound=15):
    """

    A target of 5 Hz is used, with simulations with 4, 5, and 6 spikes taking extra trials. This target and the
    range to expand the number of trials can be implemented relatively easily in the future.

    Assumptions
    - GABAa gmax per synapse = 350 pS
    - AMPA+NMDA gmax per synapse = 1000 pS [AMPA:NMDA ratio = 1]

    :param f: file handle to write results
    :type f: TextIO
    """
    step = 10
    if inh_syn_list is None:
        inh_syn_list = range(0, 800 + step, step)
    if exc_syn_list is None:
        exc_syn_list = range(0, 800 + step, step)
    if upper_bound is not None:
        assert upper_bound > 0
    else:
        upper_bound = 10000

    prev_run = f.read()

    load_file(synapse_type)

    h.tstop = 1000
    h.useCV()

    # without a 2nd load NEURON seems to complain about synapseFile being a non-variable...
    load_file(synapse_type)

    for inh_syn in inh_syn_list:
        above_bound = False
        for exc_syn in exc_syn_list:
            line = "exc_syn: {} ({})\t inh_syn: {} ({})\n".format(exc_syn, exc_weight, inh_syn, inh_weight)
            if line in prev_run:
                continue
            f.write(line)
            logger.info(line)
            if above_bound:
                # skip values of excitation when firing rate is already above
                # 15 Hz for lower values of excitation
                f.write("\n{} +- {} (n={})\n".format(16, 0, 0))
                logger.info("{} +- {} (n={})\n".format(16, 0, 0))
            else:
                # newSynapses from methods.hoc
                h.newSynapses(exc_syn, inh_syn)
                # Frequency, Noise, Weight (# synapses)
                h.inPy(inh_freq, inh_noise, inh_weight)
                h.ex(exc_freq, exc_noise, exc_weight)
                x = []
                f.write("spikes:\n")
                num_zero = 0
                run_list = list(range(num_runs))
                for run in run_list:
                    h.run()
                    x.append(h.apc.n)
                    f.write("\t {}".format(h.apc.n))
                    logger.info("\t spikes {}".format(h.apc.n))
                    if h.apc.n == 0:
                        num_zero += 1
                        if num_zero >= num_runs/2:
                            logger.info("break from too many 0's")
                            f.write("break from too many 0's")
                            break
                    elif (4 <= h.apc.n <= 6) and len(run_list) < max_num_runs:
                        # get more accurate results for values close to 5
                        run_list.append(len(run_list))
                    elif h.apc.n > upper_bound:
                        logger.info("break from above {}".format(upper_bound))
                        f.write("break from above {}".format(upper_bound))
                        above_bound = True
                        break
                mean = np.mean(x)
                std = np.std(x)
                f.write("\n{} +- {} (n={})\n".format(mean, std, run + 1))
                logger.info("{} +- {} (n={})\n".format(mean, std, run + 1))


def check_synapse_strengths(f, syn_list=(5, 100), target=5.,
                            num_runs=5, max_num_runs=10):
    """
    Find the number of synapses that yields the target (Hz)

    Assumptions:
    - Same number of synapses
    :param f: file handle to write results
    :type f: TextIO
    """

    h.tstop = 1000
    inh_freq, inh_noise = 5, 1
    exc_freq, exc_noise = 5, 1
    exc_weights = [1, 2, 5, 7, 10, 15, 20]
    inh_weights = [1, 2, 5, 7, 10, 15, 20]

    means = []
    stds = []
    weights = []

    h.inPy(1, inh_noise, 0)  # Frequency, Noise, Weight (# synapses)
    h.ex(1, exc_noise, 0)  # Frequency, Noise, Weight (uS)
    h.run()
    h.useCV()
    for syn in syn_list:
        # newSynapses from methods.hoc
        h.newSynapses(syn, syn)
        for inh_weight in inh_weights:
            above_max_f = False
            for exc_weight in exc_weights:
                f.write(
                        "exc_weight: {} ({})\t inh_weight: {} ({})\n".format(exc_weight, syn, inh_weight, syn))
                logger.info("exc_weight:{} ({})\t inh_weight: {} ({})\n".format(exc_weight, syn, inh_weight, syn))
                if above_max_f:
                    # skip values of excitation when firing rate is already above
                    # 11 Hz for lower values of excitation
                    f.write("\n{} +- {} (n={})\n".format(12, 0, 0))
                    logger.info("{} +- {} (n={})\n".format(12, 0, 0))
                else:

                    h.inPy(inh_freq, inh_noise, inh_weight)
                    h.ex(exc_freq, exc_noise, exc_weight)
                    x = []
                    f.write("spikes:\n")
                    num_zero = 0
                    run_list = list(range(num_runs))
                    for run in run_list:
                        h.run()
                        x.append(h.apc.n)
                        f.write("\t {}".format(h.apc.n))
                        logger.info("\t spikes {}".format(h.apc.n))
                        if not (0 < h.apc.n <= 11):
                            if h.apc.n == 0:
                                num_zero += 1
                                if num_zero >= num_runs/5:
                                    logger.info("break from too many 0's")
                                    f.write("break from too many 0's")
                                    break
                            else:
                                logger.info("break from above 11")
                                f.write("break from above 11")
                                above_max_f = True
                                break
                        elif (4 <= h.apc.n <= 6) and len(run_list) < max_num_runs:
                            # get more accurate results for values close to 5
                            run_list.append(len(run_list))
                    mean = np.mean(x)
                    std = np.std(x)
                    f.write("\n{} +- {} (n={})\n".format(mean, std, run + 1))
                    logger.info("{} +- {} (n={})\n".format(mean, std, run + 1))
                    means.append(mean)
                    stds.append(std)
                    weights.append((exc_weight, inh_weight))

    arr = np.abs(target - np.array(means))
    idx = np.argmin(arr)
    return weights[idx]
