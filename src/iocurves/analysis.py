import ast
import os

import logging
import numpy as np
import pandas as pd


logger = logging.getLogger("iocurves analysis")


def boltzman(x, xmid, tau):
    """
    evaluate the boltzman function with midpoint xmid and time constant tau over x
    """
    return 1./(1. + np.exp(-(x - xmid)/tau))


def sigmoid(x, x0, k):
    """
    evaluate sigmoid function slope k and midpoint x0 over x
    """
    y = 1/(1 + np.exp(-k*(x - x0)))
    return y


# keep some results in memory for quick re-evaluation (ast can take a while)
memoised = dict()


def get_params(save_name):
    """Convert a name to python data base don a known structure.

    :param save_name: str in format:
        `"[file name[_KCC2]]_[synapse type]_[synapse numbers]_[syn input]_...[recording location]_[trials].txt"`
    :type save_name: str
    :return: file_name, synapse_type, synapse_numbers, syn_input, non_default_dict, location, trials
    :rtype: (str, int, list of int, dict of str, dict, str, int)
    """
    if save_name in memoised:
        return memoised[save_name]
    params = save_name.split('_')
    if params[1] == 'KCC2':
        params[0] = params[0] + '_' + params[1]
        del params[1]
    file_name, synapse_type, synapse_numbers, syn_input = params[:4]
    location, trials = params[-2:]
    non_default_dict = {"diam": None, "pa_kcc2": None}
    if len(params) > 6:
        extra_keys = list(non_default_dict.keys())
        for i in range(4, len(params) - 2):
            non_default_dict[extra_keys[i - 4]] = params[i]
    # process some of the params
    synapse_numbers = ast.literal_eval(synapse_numbers)
    syn_input = ast.literal_eval(syn_input)
    for key, value in non_default_dict.items():
        try:
            non_default_dict[key] = ast.literal_eval(non_default_dict[key])
        except ValueError:
            # do nothing
            pass
    memoised[save_name] = (file_name, synapse_type, synapse_numbers, syn_input, non_default_dict, location, trials)
    return file_name, synapse_type, synapse_numbers, syn_input, non_default_dict, location, trials


def get_var(df, recorded_var, mean=True):
    """
    Get variable from dataframe one level deep and optionally the mean.

    :param df: Recorded data over time
    :type df: pd.DataFrame
    :param recorded_var: Variable to retrieve from df
    :type recorded_var: str
    :param mean: Include the mean in the return tuple
    :type mean: bool
    :return: Dataframe of variable and Series of the mean for the variable for each time step
    :rtype: (pd.DataFrame, pd.Series or None)
    """
    all_var_columns = df.xs(recorded_var, level=1, axis=1)
    if mean:
        # mean over the columns
        mean_val = all_var_columns.mean(axis=1)
    else:
        mean_val = None
    return all_var_columns, mean_val


def moving_average(spike_indices, time_bin=1., backward=True):
    """ Calculate the moving average filter for the spike_indices according to time_bin (in seconds).

    :param spike_indices: Array of spike/no spike boolean type over time.
    :type spike_indices: np.ndarray of bool
    :param time_bin: Size of sliding window computation (in seconds)
    :type time_bin: float
    :param backward: Return values are for t-time_bin (True) or t+time_bin (False)
    :type backward: bool
    :return: Instantaneous firing rate for each point in time (overlapping windows)
    :rtype: np.ndarray or float
    """
    from src.iocurves.sim import DT
    time_bin_size = int(time_bin*1000/DT)
    ifr = np.cumsum(spike_indices, dtype=float)
    if backward:
        ifr[:-time_bin_size] = ifr[time_bin_size:] - ifr[:-time_bin_size]  # backward window
        # amend end (if backward window) of IFR trace to be flat
        ifr[-time_bin_size:] = ifr[-time_bin_size - 1]
    else:
        ifr[time_bin_size:] = ifr[time_bin_size:] - ifr[:-time_bin_size]  # forward window
        ifr = ifr[time_bin_size - 1:]
    return ifr/time_bin


def get_inst_firing_rate(spike_arr, time_bin=1., slide=True):
    """

    :param spike_arr:
    :type spike_arr: np.ndarray of int or list of int
    :param time_bin: Size of sliding window computation (in seconds)
    :type time_bin: float
    :param slide:
    :type slide: bool
    :return:
    :rtype:
    """
    from src.iocurves.sim import DT

    if type(spike_arr) != np.ndarray:
        spike_arr = np.array(spike_arr)
    spike_indices = np.diff(spike_arr) > 0
    spike_indices = np.append(spike_indices, [False])  # add an element due to diff losing one
    if slide:
        # super fast
        ifr = moving_average(spike_indices, time_bin)
    else:
        ifr = np.zeros(shape=len(spike_indices))
        time_bin_size = int(time_bin*1000/DT)
        # slower for moving average, but fast enough for non-overlapping intervals
        for i in range(0, int(len(spike_indices) - time_bin_size), time_bin_size):
            ifr[i:i + time_bin_size] = sum(spike_indices[i:i + time_bin_size])/time_bin

    return ifr


def save_to_file(title, result):
    path = os.path.join('results', title)
    if type(result) is list:
        np.save(path + '.npy', result)
    elif type(result) is pd.DataFrame:
        result.to_hdf(path, 'table')
    logger.info("saved")


def load_from_file(title):
    """

    :param title:
    :return: list() or None
    """
    path = os.path.join('results', title)
    try:
        return pd.read_hdf(path, 'table')
    except IOError:
        try:
            return np.load(path+".npy", allow_pickle=True)
        except OSError:
            return None


def get_data(cl_state_trials, ifr_windowsize, time_points, var='spikes'):
    from src.iocurves.sim import DT
    FRdf = pd.DataFrame(index=time_points)
    exc = set()
    inh = set()
    for key, df in cl_state_trials.items():
        file_name, synapse_type, synapse_numbers, syn_input, diam, location, trials = get_params(key)
        spike_results, spike_mean = get_var(df, var)
        # number_spikes = spike_mean.iloc[-1]
        # convert spike times to instantaneous firing rate
        trial_length = spike_results.shape[1] + 1
        if var == 'spikes':
            for j in range(1, trial_length):
                spike_results.loc[:, j] = get_inst_firing_rate(spike_results[j], time_bin=ifr_windowsize)
        spike_mean.loc[:] = np.mean(spike_results, axis=1)
        if time_points is None:
            # no specific time points, use all time points
            ifr = spike_mean
            trial = pd.DataFrame({
                (syn_input['in'], syn_input['ex']): ifr
                }, index=spike_mean.index)
        else:
            ifr = []
            for i, time_point in enumerate(time_points):
                start = int(time_point/DT)
                end = int(start + ifr_windowsize*1000/DT)
                ifr.append(np.mean(spike_mean.iloc[start:end]))
                # logger.info(ifr[i])
            trial = pd.DataFrame({
                (syn_input['in'], syn_input['ex']): ifr
                }, index=time_points)
        FRdf = pd.concat([FRdf, trial], axis=1)
        exc.add(syn_input['ex'])
        inh.add(syn_input['in'])
    return FRdf, exc, inh
