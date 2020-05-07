import itertools

import logging
import numpy as np

from src.iocurves.analysis import get_var
from src.iocurves.sim import do_runs

logger = logging.getLogger("run range")


def run_precise(exc_input=None, inh_input=None, *args, **kwargs):
    # type: (list, list, list, dict) -> (dict,dict)
    if exc_input is None:
        exc_input = [0, 0.01, 0.05, 0.1, 0.5, 1, 2, 4, 10, 20, 40, 50, 100, 200, 500, 1000]
    if inh_input is None:
        inh_input = list(exc_input)
    syn_nums = list(itertools.product(inh_input, exc_input))
    kcc2_runs = {}
    static_runs = {}

    for i, syn_nums_pair in enumerate(syn_nums):
        logger.info("="*30)
        logger.info(syn_nums_pair)
        result, save_name, result_kcc2, save_name_kcc2 = do_runs(*args,
                                                                 syn_input={
                                                                     'ex': syn_nums_pair[1],
                                                                     'in': syn_nums_pair[0]
                                                                     },
                                                                 **kwargs)
        kcc2_runs[save_name_kcc2] = result_kcc2
        static_runs[save_name] = result
    return static_runs, kcc2_runs


def run_dynamic(lower_ex, upper_ex, lower_inh, upper_inh, mid='log', fr_diff=10, syn_diff=1.0, **common_args):
    if lower_ex <= 0 and mid == 'log':
        lower_ex = 1
        logger.info('lower_ex change to {}'.format(lower_ex))

    def run_initial(lower_ex, upper_ex, lower_inh, upper_inh, **common_args):
        kcc2_runs_initial = {}
        static_runs_initial = {}
        result_leli, save_name_leli, result_kcc2_leli, save_name_kcc2_leli = do_runs(syn_input={
            'ex': lower_ex,
            'in': lower_inh
            },
                **common_args)
        result_ueli, save_name_ueli, result_kcc2_ueli, save_name_kcc2_ueli = do_runs(syn_input={
            'ex': upper_ex,
            'in': lower_inh
            },
                **common_args)
        result_leui, save_name_leui, result_kcc2_leui, save_name_kcc2_leui = do_runs(syn_input={
            'ex': lower_ex,
            'in': upper_inh
            },
                **common_args)
        result_ueui, save_name_ueui, result_kcc2_ueui, save_name_kcc2_ueui = do_runs(syn_input={
            'ex': upper_ex,
            'in': upper_inh
            },
                **common_args)
        spikes_leli, spikes_mean_leli = get_var(result_leli, 'spikes', True)
        spikes_ueli, spikes_mean_ueli = get_var(result_ueli, 'spikes', True)
        spikes_leui, spikes_mean_leui = get_var(result_leui, 'spikes', True)
        spikes_ueui, spikes_mean_ueui = get_var(result_ueui, 'spikes', True)
        lower_spikes = max(spikes_mean_leli)
        upper_spikes = max(spikes_mean_ueli)

        kcc2_runs_initial[save_name_kcc2_leli] = result_kcc2_leli
        kcc2_runs_initial[save_name_kcc2_ueli] = result_kcc2_ueli
        kcc2_runs_initial[save_name_kcc2_leui] = result_kcc2_leui
        kcc2_runs_initial[save_name_kcc2_ueui] = result_kcc2_ueui
        static_runs_initial[save_name_leli] = result_leli
        static_runs_initial[save_name_ueli] = result_ueli
        static_runs_initial[save_name_leui] = result_leui
        static_runs_initial[save_name_ueui] = result_ueui

        return static_runs_initial, kcc2_runs_initial, lower_spikes, upper_spikes

    static_runs, kcc2_runs, lower_spikes, upper_spikes = run_initial(lower_ex, upper_ex, lower_inh, upper_inh,
                                                                     **common_args)
    syn_nums = [(lower_ex, upper_ex, lower_inh, upper_inh, lower_spikes, upper_spikes, True)]
    seen = set()
    for syn_nums_pair in syn_nums:
        logger.info("="*60 + "len(seen) = {}".format(len(seen)))
        logger.info(syn_nums_pair)
        lower_ex, upper_ex, lower_inh, upper_inh, lower_spikes, upper_spikes, increase = syn_nums_pair
        if ((lower_ex, upper_ex, lower_inh, upper_inh) in seen) or (
                syn_diff is not None and (
                (upper_ex - lower_ex < syn_diff) or (upper_inh - lower_inh < syn_diff))):
            continue

        if len(seen) > 500:
            logger.info("BREAK")
            break
        seen.add((lower_ex, upper_ex, lower_inh, upper_inh))

        if mid == 'half':
            mid_ex = lower_ex + int((upper_ex - lower_ex)/2)
            mid_inh = lower_inh + int((upper_inh - lower_inh)/2)
        elif mid == 'log':
            mid_ex = 10 ** (np.log10(upper_ex/lower_ex)/2)
            mid_inh = 10 ** (np.log10(upper_inh/lower_inh)/2)

        result_meli, save_name_meli, result_kcc2_meli, save_name_kcc2_meli = do_runs(syn_input={
            'ex': mid_ex,
            'in': lower_inh
            },
                **common_args)
        result_lemi, save_name_lemi, result_kcc2_lemi, save_name_kcc2_lemi = do_runs(syn_input={
            'ex': lower_ex,
            'in': mid_inh
            },
                **common_args)
        result_memi, save_name_memi, result_kcc2_memi, save_name_kcc2_memi = do_runs(syn_input={
            'ex': mid_ex,
            'in': mid_inh
            },
                **common_args)
        result_uemi, save_name_uemi, result_kcc2_uemi, save_name_kcc2_uemi = do_runs(syn_input={
            'ex': upper_ex,
            'in': mid_inh
            },
                **common_args)
        result_meui, save_name_meui, result_kcc2_meui, save_name_kcc2_meui = do_runs(syn_input={
            'ex': mid_ex,
            'in': upper_inh
            },
                **common_args)
        spikes_meli_df, spikes_mean_meli_series = get_var(result_meli, 'spikes', True)
        spikes_lemi_df, spikes_mean_lemi_series = get_var(result_lemi, 'spikes', True)
        spikes_memi_df, spikes_mean_memi_series = get_var(result_memi, 'spikes', True)
        spikes_uemi_df, spikes_mean_uemi_series = get_var(result_uemi, 'spikes', True)
        spikes_meui_df, spikes_mean_meui_series = get_var(result_meui, 'spikes', True)

        spikes_mean_meli = max(spikes_mean_meli_series)
        spikes_mean_lemi = max(spikes_mean_lemi_series)
        spikes_mean_memi = max(spikes_mean_memi_series)
        spikes_mean_uemi = max(spikes_mean_uemi_series)
        spikes_mean_meui = max(spikes_mean_meui_series)

        if spikes_mean_meli < upper_spikes - fr_diff:
            syn_nums.append((mid_ex, upper_ex, lower_inh, upper_inh, spikes_mean_meli, upper_spikes, True))
        if spikes_mean_meli > lower_spikes + fr_diff:
            syn_nums.append((lower_ex, mid_ex, lower_inh, upper_inh, lower_spikes, spikes_mean_meli, False))

        if fr_diff < spikes_mean_lemi < upper_spikes - fr_diff:
            syn_nums.append((lower_ex, upper_ex, mid_inh, upper_inh, spikes_mean_lemi, upper_spikes, True))
        if spikes_mean_lemi > lower_spikes + fr_diff:
            syn_nums.append((lower_ex, upper_ex, lower_inh, mid_inh, lower_spikes, spikes_mean_lemi, False))

        if fr_diff < spikes_mean_uemi < upper_spikes - fr_diff:
            syn_nums.append((lower_ex, upper_ex, mid_inh, upper_inh, spikes_mean_uemi, upper_spikes, True))
        if spikes_mean_uemi > lower_spikes + fr_diff:
            syn_nums.append((lower_ex, upper_ex, lower_inh, mid_inh, lower_spikes, spikes_mean_uemi, False))

        if increase and spikes_mean_meui < upper_spikes - fr_diff:
            syn_nums.append((mid_ex, upper_ex, lower_inh, upper_inh, spikes_mean_meui, upper_spikes, True))
        if spikes_mean_meui > lower_spikes + fr_diff:
            syn_nums.append((lower_ex, mid_ex, lower_inh, upper_inh, lower_spikes, spikes_mean_meui, False))

        kcc2_runs[save_name_kcc2_meli] = result_kcc2_meli
        kcc2_runs[save_name_kcc2_lemi] = result_kcc2_lemi
        kcc2_runs[save_name_kcc2_memi] = result_kcc2_memi
        kcc2_runs[save_name_kcc2_uemi] = result_kcc2_uemi
        kcc2_runs[save_name_kcc2_meui] = result_kcc2_meui
        static_runs[save_name_meli] = result_meli
        static_runs[save_name_lemi] = result_lemi
        static_runs[save_name_memi] = result_memi
        static_runs[save_name_uemi] = result_uemi
        static_runs[save_name_meui] = result_meui
    return static_runs, kcc2_runs
