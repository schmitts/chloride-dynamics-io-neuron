import logging

import pandas as pd
from neuron import h
from src.iocurves.analysis import load_from_file, save_to_file
from src.utils.nrnpy import get_base_vm_cli, load_file
from src.utils.timing import current_time

logger = logging.getLogger("cl time")

TSTOP = 10000.0  # ms
SUBSAMPLE = 200  # don't load every point into memory
DT = 0.025*SUBSAMPLE


def pyrun(file_name, synapse_type=1, synapse_numbers=(100, 100), syn_input=None, diam=None, pa_kcc2=None,
          location='axon', trials=1, save=True, tstop=TSTOP, **kwargs):
    """
    Run a NEURON simulation for a neuron specified in ``file_name`` with input specified by other parameters provided.

    :param file_name: Neuron definition (excluding '.hoc').
    :type file_name: str
    :param synapse_type: Type of synapses to use (0 for frequency-based 'f-in', 1 for persistent conductance, 'gclamp').
    :type synapse_type: int
    :param synapse_numbers: Number of (E, I) on the neuron.
    :type synapse_numbers: (int, int)
    :param syn_input: Mapping of excitatory/inhibitory type to input strength. {'ex', E, 'in: I}
    :type syn_input: dict[str, int]
    :param diam: Re-specify diam for specific regions of the neuron. Valid: 'ldend', 'bdend', 'soma', 'axon'.
    :type diam: dict[str: float]
    :param pa_kcc2: Strength of KCC2.
    :type pa_kcc2: float
    :param location: Location to recording firing rate.
    :type location: str
    :param trials: Number of repeated simulations to run.
    :type trials: int
    :param save: Whether to load/save the results from/to file.
    :type save: bool
    :param tstop: Length of simulation (ms).
    :type tstop: float
    :param kwargs: Other keywords are ignored.

    :return: Pair of DataFrame with results and name of save file (even if not saved, the name is generated).
    :rtype: (pd.DataFrame, str)
    """
    if syn_input is None:
        syn_input = {'in': 5, 'ex': 5}
    save_name = "{}_{}_{}_{}".format(file_name, synapse_type, synapse_numbers, syn_input)
    save_name += "_{}".format(diam) if diam is not None else ''
    save_name += "_{}".format(pa_kcc2) if pa_kcc2 is not None else ''
    save_name += "_{}_{}".format(location, trials)
    logger.info(save_name)
    if save:
        loaded = load_from_file(save_name)
        if loaded is not None:
            return loaded, save_name

    load_file(file_name)
    load_file(file_name)
    if diam is not None:
        for seg in diam.keys():
            nrn_seg = get_compartment(seg)
            nrn_seg.diam = diam[seg]
    if pa_kcc2 is not None:
        hoc_cmd = "forall {" + """
            Pa_KCC2 = {pa_kcc2}e-5""".format(pa_kcc2=pa_kcc2) +\
                  " }"
        h(hoc_cmd)
    compartment = get_compartment(location)
    if file_name.find('distal') > -1:
        cli_rec_loc = get_compartment('ldend')
    elif file_name.find('proximal') > -1 or file_name.find('proximal') > -1:
        cli_rec_loc = get_compartment('bdend')
    else:
        cli_rec_loc = get_compartment('soma')

    logger.info("recording cli from {}".format(cli_rec_loc.hname()))

    h("access {}".format(location))
    h.changeSynapseType(synapse_type)
    h.newSynapses(synapse_numbers[0], synapse_numbers[1])
    h.inPy(0)
    h.ex(0)
    vm_init, cli = get_base_vm_cli(file_name, compartment)
    h.v_init = vm_init
    h_str = """
        forall{""" + """
            cli = {0}
            cli0_cl_ion = {0}
            """.format(cli)
    if file_name.find("KCC2") > 0:
        h_str += """cli0_KCC2 = {0}
        """.format(cli)
    h_str += "}"
    h(h_str)

    h.tstop = tstop

    if synapse_type == 0:
        # Hz, duration (s), start (ms), noise, weight/channels
        h.inPy(syn_input['in'], h.tstop/1000, 0, 1, 1)
        h.ex(syn_input['ex'], h.tstop/1000, 0, 1, 1)
    else:
        h.inPy(syn_input['in'])
        h.ex(syn_input['ex'])

    # create recording vector objects
    t_rec = h.Vector()
    v_rec = h.Vector()
    cli_rec = h.Vector()
    spike_rec = h.Vector()

    t_rec.record(h._ref_t)
    v_rec.record(compartment(0.5)._ref_v)
    cli_rec.record(cli_rec_loc(0.5)._ref_cli)
    spike_rec.record(h.apc._ref_n)

    time_past = current_time('ms')
    logger.info("using {}...".format(file_name))
    assert trials > 0
    logger.info(save_name)
    logger.info("trial # | # spikes")
    df = pd.DataFrame()
    for i in range(trials):
        trial_num = i + 1
        h.run()
        logger.info("{:7} | {:8}".format(trial_num, h.apc.n))

        temp_dict = {
            (trial_num, 'v'):      v_rec.as_numpy(),
            (trial_num, 'cli'):    cli_rec.as_numpy(),
            (trial_num, 'spikes'): spike_rec.as_numpy()
            }
        recording = pd.DataFrame(temp_dict, index=t_rec.to_python())

        df = pd.concat([df, recording], axis=1)

    logger.info("time taken: {}ms".format(current_time('ms') - time_past))

    if save:
        save_to_file(save_name, df)
    return df, save_name


GC = 0


def do_runs(file_name, plot=[()] or True, syn_input=None, ifr_windowsize=1., **kwargs):
    """
    Run a neuron with and without KCC2 as specified in ``file_name`` with synaptic input as in ``syn_input``.

    Specific inputs can be plotted by specifying a list of (E,I) inputs to plot.


    :param file_name: Neuron definition (excluding '.hoc').
    :type file_name: str
    :param plot: Whether to plot static vs dynamic Cl-.
        If a list of tuples, will only plot if (syn_input['ex'], syn_input['in']) is present.
    :type plot: bool or List[Tuple[float, float]]]
    :param syn_input: Excitatory "ex" and inhibitory "in" input, specified as a mapping.
        E.g. {'ex': 5, 'in': 5}
    :type syn_input: dict[str, float]
    :param ifr_windowsize: Instantaneous firing rate window size. Used for plotting, ignored otherwise.
    :type ifr_windowsize: float
    :param kwargs: Other keywords to pass to :meth:`pyrun`.
    :return:
    :rtype:
    """
    global GC
    result, save_name = pyrun(file_name, syn_input=syn_input, **kwargs)
    result_kcc2, save_name_kcc2 = pyrun(file_name + "_KCC2", syn_input=syn_input, **kwargs)
    if plot or (type(plot) is dict and plot.count((syn_input['ex'], syn_input['in'])) > 0):
        from src.iocurves.vis import plot_compare_dcl
        plot_compare_dcl(save_name, [result, result_kcc2], ifr_windowsize=ifr_windowsize, **kwargs)

    result = result.iloc[::SUBSAMPLE]
    result_kcc2 = result_kcc2.iloc[::SUBSAMPLE]

    # try clear some RAM
    if GC == 10:
        import gc
        gc.collect()
        GC = 0
    else:
        GC += 1

    return result, save_name, result_kcc2, save_name_kcc2


is_menu = 0


def show_menu():
    global is_menu
    if is_menu == 0:
        h.nrnmainmenu()  # create main menu
        h.nrncontrolmenu()  # crate control menu
        is_menu = 1


def get_compartment(recording_location='soma'):
    if recording_location == 'soma':
        compartment = h.soma
    elif recording_location == 'axon':
        compartment = h.axon
    elif recording_location == 'ldend':
        compartment = h.ldend
    elif recording_location == 'bdend':
        compartment = h.bdend
    else:
        raise Exception("specify valid location")
    return compartment
