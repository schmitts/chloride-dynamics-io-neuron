# coding=utf-8
"""
CONFIGURATION TYPES
"""
from neuron import h

from src.utils.nrnpy import load_file


def load_config():
    """Base configuration method."""
    load_file("config.hoc", set_file_name=False, force_reload=False)


def nruns(num_runs=1):
    """Set number of repeated trials."""
    load_config()
    h("nRepeats_config = {}".format(num_runs))


def balanced_input(index_synced=1):
    """For balanced input, only run simulations where exc input = inh input."""
    load_config()
    h("""index_synced={}""".format(index_synced))


def config_synapses(exc, inh):
    """Set frequency of synapses."""
    load_config()
    h.nexc_hz_config = exc
    h.ninh_hz_config = inh
    # h.set_optimal_weights(exc, inh)


def single_hz_and_control(hz=5, control=True):
    """Overwrite frequency list in config_stims.hoc so only a single frequency is in the arrays.
    Useful when increasing the number of synapses while keeping input frequency the same."""
    load_config()
    single_stim_freq_proc = """
        proc defineSingleStimFreq(){
            singleStimFreq = new Vector()
        """
    if control:
        single_stim_freq_proc +=\
            "singleStimFreq.append(0)"
    single_stim_freq_proc += """
            singleStimFreq.append({0})
                        """.format(hz)
    single_stim_freq_proc += "}"
    h(single_stim_freq_proc)
    h("""
        proc proximalHzConfig(){
            print "proximalHzConfig custom Python"
            proximalWeights()
            defineSingleStimFreq()
            setInhibitStims(singleStimFreq)
            setExciteStims(singleStimFreq)
        }
        proc distalHzConfig(){
            print "distalHzConfig custom Python"
            distalWeights()
            defineSingleStimFreq()
            setInhibitStims(singleStimFreq)
            setExciteStims(singleStimFreq)
        }
    """)
    print("single_hz_and_control applied")


def config_persistent_synapses(**kwargs):
    """Force synapses to be fluctuating conductances."""
    load_config()
    h("""
        persistentSynapses_config = 1 // 0 for none, 1 for both, 2 for inh, 3 for exc
        proc proximalInhPersistantConfig(){
            print "proximalInhPersistantConfig"
            proximalInhPersistant()
        }
        proc proximalExcPersistantConfig(){
            print "proximalExcPersistantConfig"
            excPersistant()
        }
        proc distalInhPersistantConfig(){
            print "distalInhPersistantConfig"
            distalLowerInhPersistant()
        }
        proc distalExcPersistantConfig(){
            print "distalExcPersistantConfig"
            excPersistant()
        }
    """)
    # other hoc assignments
    for kw, val in kwargs.items():
        h("{} = {}".format(kw, val))
