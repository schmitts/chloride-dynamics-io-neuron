# coding=utf-8
"""
Comparisons to set up and run in NEURON using `hoc_files`

Each comparison is called via the `nrnpy.start` method.
Thus, comparisons require arguments for `filename` and `param_list` (wthe structure of which can be different for
each method)

"""
import numpy as np
from neuron import h

from src.utils.nrnconfig import load_config
from src.utils.nrnpy import get_base_vm_cli, load_file, save_run

SKIP_IF_EXISTS = True


# ----------------------------------------------------------------------------------------------------------------------
# TEST TYPES
# ----------------------------------------------------------------------------------------------------------------------

def just_run(filename, param_list=None):
    """
    Proxy for save_run, that accepts a param_list as expected in `nrnpy.start`
    """
    load_file(filename)
    if param_list is not None:
        if type(param_list) is str:
            h(param_list)
        elif type(param_list[0]) is str:
            for cmd in param_list:
                h(cmd)
    save_run(filename)


def compare_pCl(filename, param_list=None):
    """Change leak chloride permeability (ratio to leak potassium conductance)"""
    pcl_list = param_list
    if not pcl_list:
        pcl_list = [0.1, 0.4, 0.8, 1.2]
    print("pcl_list:{}".format(pcl_list))
    load_file(filename)
    # list of (exc,inh) synapses numbers
    for pcl in pcl_list:
        h.changePCl(pcl)
        save_run("{}_{}".format(filename, pcl))


def compare_pas(filename, param_list=None):
    """Change input resistance by adjusting the passive leak channel.
    The leak channel has K+, Na+, and Cl- conductances with a set ratio (see config.hoc).
    Calling h.changePas(<new passive K conductance>) respects the ratio.
    Saved filenames include the relative change in PAS, the Input Resistance, the steady-state cli (at 0 input),
    and membrane potential.
    """
    rel_gkpbar = param_list
    if not rel_gkpbar:
        rel_gkpbar = [0.8, 1, 1.2, 1.5]
    print("rel_gkpbar:{}".format(rel_gkpbar))
    load_file(filename)
    # list of (exc,inh) synapses numbers
    imp = h.Impedance()
    imp.loc(0.5, sec=h.soma)
    for p in rel_gkpbar:
        h.changePas(h.g_pas_k_config*p)
        vm, cli = get_base_vm_cli(filename, compartment=h.soma)
        imp.compute(0, 1)  # calculate impedance at 0 Hz (i.e. resistance)
        ir = imp.input(0.5, sec=h.soma)  # input resistance at 0.5
        print("when p={}, Rm={}, [Cl]i={}, and Vm={}".format(p, ir, cli, vm))
        save_run("{}_[{},{:.2f},{:.2f},{:.2f}]".format(filename, p, ir, cli, vm))


def compare_pkcc2(filename, param_list=None):
    """Change pump-strength of KCC2"""
    pkcc2_list = param_list
    if not pkcc2_list:
        base = 1.9297e-5  # default value from mod file
        pkcc2_list = [base/2, base, base*2, base*4]
    print("pkcc2_list:{}".format(pkcc2_list))
    load_file(filename)
    for pkcc2 in pkcc2_list:
        cmd = "forall{" +\
              "Pa_KCC2 = {}".format(pkcc2) +\
              "}"
        h(cmd)
        save_run("{}_{}".format(filename, pkcc2))


def compare_pkcc2_homo(filename, param_list=None):
    """Change homogeneity of KCC2 pump strength"""
    axon_pkcc2_list = [1]
    ldend_pkcc2_list = [0.5, 1, 2, 4]
    if param_list is not None:
        if np.iterable(param_list[0]) and type(param_list[0]) is not str:
            ldend_pkcc2_list, axon_pkcc2_list = param_list
        else:
            ldend_pkcc2_list = param_list
    base = 1.9297e-5
    print("ldend_pkcc2_list:{}".format(ldend_pkcc2_list))
    load_file(filename)

    if "KCC2" not in filename:
        # run a single simulation without KCC2.
        savename = "{}_({:.2f},{:.2f})".format(filename, 0, 0)
        save_run(savename)
        return

    cmd_bdend = "bdend.Pa_KCC2 = {}".format(base)
    h(cmd_bdend)
    for pkcc2 in ldend_pkcc2_list:
        cmd = "ldend.Pa_KCC2 = {}".format(base*pkcc2)
        print(cmd)
        h(cmd)
        for axon_pkcc2 in axon_pkcc2_list:
            axon_cmd = "axon.Pa_KCC2 = {}".format(base*axon_pkcc2)
            print(axon_cmd)
            h(axon_cmd)
            save_run("{}_({:.2f},{:.2f})".format(filename, pkcc2, axon_pkcc2))


def compare_diam(filename, param_list=None):
    """Change distal dendrite's diameter"""
    diam_list = param_list
    if not diam_list:
        diam_list = [0.5, 1, 1.5, 2]
    print("diam_list:{}".format(diam_list))
    load_file(filename)
    for diam in diam_list:
        h("""ldend.diam = {}""".format(diam))
        save_run("{}_{}".format(filename, diam))


def compare_cli(filename, param_list=None):
    """Change internal chloride concentration (for constant/static) runs"""
    cli_list = param_list
    if not cli_list:
        cli_list = [2.5, 5, 7.5, 10, 12.5, 15, 20, 25, 30]
    print("cli_list:{}".format(cli_list))
    load_file(filename)
    for cli in cli_list:
        # chloride is set inside run.hoc through save_run --> doRun(set_cli)
        cmd = "forall{" + """
                          cli={cli}
                          cli0_cl_ion={cli}
                          cli0_KCC2={cli}
                       """.format(cli=cli) +\
              "}"
        h(cmd)
        save_run("{}_{}".format(filename, cli), set_cli=cli)


def compare_duration(filename, param_list=None):
    """Change tstop to compare"""
    tstop_list = param_list
    if not tstop_list:
        base = h.tstop_config
        tstop_list = [base/10, base/2, base, base*2, base*10]
    print("tstop_list:{}".format(tstop_list))
    load_file(filename)
    for tstop in tstop_list:
        h.tstop_config = tstop  # doRun() calles loadParameters() in run.hoc which retrieves config variables
        h.tstop = tstop  # we set this here anyway
        save_run("{}_{}".format(filename, tstop))


def compare_dynamic_K(filename, param_list=None):
    """Run simulations with dynamic potassium (True) or not (False)
    The base filename should be used (e.g. 'distal') such that
    - distal is the baseline
    - distal_KCC2_True is KCC2 WITH dynamic potassium and
    - distal_KCC2_False is KCC2 WITHOUT dynamic potassium
    """
    if 'KCC2' in filename:
        print("ignoring {} in 'compare_dynamic_K' because param_list should cover this case".format(filename))
        return
    K_list = param_list
    if not K_list:
        K_list = ['control', True, False]
    print("K_list:{}".format(K_list))

    for dyn_K in K_list:
        fname = filename if dyn_K == 'control' else filename + "_KCC2"
        assert (dyn_K != 'control' and 'KCC2' in fname) or (dyn_K == 'control' and 'KCC2' not in fname)
        load_file(fname)
        # get steady state here as it isn't done in run.hoc anymore with the set_cli arg
        vm, cli = get_base_vm_cli(fname, compartment=h.ldend)
        h.v_init = vm
        if dyn_K == 'control':
            print("control (no KCC2 nor KCC2_pot")
            cmd = "forall{" + """
                              cli={cli}
                              cli0_cl_ion={cli}
                           """.format(cli=cli) +\
                  "}"
            h(cmd)
            save_run("{}_{}".format(fname, 0), set_cli=cli)
            pot = ""
        elif dyn_K:
            # then add KCC2 potassium
            cmd = """forall{
                uninsert pasghk
                uninsert KCC2
                insert pas
                insert KCC2_pot
                g_pas = 0.00021
                }"""
            h(cmd)
            pot = "_pot"
        else:
            # already have KCC2 in the neuron
            pot = ""
        cmd = "forall{" + """
                  cli={cli}
                  cli0_cl_ion={cli}
                  cli0_KCC2{pot}={cli}
               """.format(cli=cli, pot=pot) +\
              "}"
        h(cmd)
        save_run("{}_{}".format(fname, int(dyn_K)), set_cli=cli)


def compare_synapses(filename, param_list=None, errors=0):
    """Run simulations with different number of (E,I) synapses"""
    if type(param_list) is tuple:
        skip, syn_list = param_list
    else:
        syn_list = param_list
        skip = SKIP_IF_EXISTS
    if not syn_list:
        syn_list = [(120, 10), (140, 50), (180, 130)]
    print("syn_list:{}".format(syn_list))
    load_file(filename)
    vm, cli = get_base_vm_cli(filename)
    h.v_init = vm
    for i, (exc, inh) in enumerate(syn_list):
        h.newSynapses(exc, inh)
        try:
            save_run("{}_{},{}".format(filename, exc, inh), set_cli=cli, skip_if_exists=skip)
        except RuntimeError as re:
            print("Error {} occured on E={} & I={}".format(re, exc, inh))
            print("leftover params = {}".format(syn_list[i:]))
            if errors == 2:
                print("second time error occured, bailing")
                import sys
                sys.exit(-1)
            else:
                compare_synapses(filename, param_list=syn_list[i:], errors=errors + 1)


def compare_weights(filename, param_list=None):
    """Run Simulations with different strengths of synapses (keeping the actual number from load_config the same.
    """
    if param_list:
        # keep param_list a list, but make first (only) arg a dict
        weights_dict = param_list[0]
        assert weights_dict['inh_weights']
        assert weights_dict['exc_weights']
    else:
        weights_dict = {
            'exc_weights': [10],
            'inh_weights': [10]
            }
    print("weight_dict:{}".format(weights_dict))
    load_file(filename)
    # load config
    load_config()
    # define hoc-accessible Vectors (lists)
    h("objref inhibWeights_config,exciteWeights_config")
    # convert python lists to NEURON Vectors
    h.inhibWeights_config = h.Vector(weights_dict['inh_weights'])
    h.exciteWeights_config = h.Vector(weights_dict['exc_weights'])
    # assign weight vectors during configWeights (which is called from (if not persistent synapses)
    # run.hoc > doRun() > loadParameters() > protocolStimVectors() > [proximal|distal]HzConfig > configWeights()
    h("""
        proc configWeights(){
            print "configWeights custom Python"
            inhibWeights = inhibWeights_config
            exciteWeights = exciteWeights_config
        }
    """)

    # iterating over weights is already done in run.hoc
    save_run("{}".format(filename))


def compare_synapses_clumped(filename, param_list=None):
    """Run simulations with synapses located close together. Clumping implemented in `methods.hoc`
    """
    syn_list = param_list
    if not syn_list:
        syn_list = [(120, 10), (140, 50), (180, 130)]
    print("syn_list:{}".format(syn_list))
    load_file(filename)
    for (exc, inh, clump) in syn_list:
        h.newSynapses(exc, inh)
        h.clumpSynapses(clump)
        save_run("{}_{},{}".format(filename, exc, inh))


def compare_synapse_types(filename, param_list=None):
    """Run simulations with different synapse types. Note that the strength of these synapses will be wildy off.

    * 0 - frequency-based ('f-in') synapses (receives input via NETSTIM)
    * 1 - fluctuating conductance 'gclamp' synapses
    * 2 - exc is 'f-in', inh is 'gclamp'
    * 3 - inh is 'f-in', exc is 'gclamp'
    """
    syn_type_list = param_list
    if not syn_type_list:
        syn_type_list = [0, 1, 2, 3]
    print("syn_type_list:{}".format(syn_type_list))
    load_file(filename)
    for type_num in syn_type_list:
        h.changeSynapseType(type_num)
        save_run("{}_{}".format(filename, type_num))


def synapses_vs_weights(filename="distal", param_list=None):
    """Clump synapses and vary synapse numbers and synapse weights.

    param_list contains tuples of (exc_syn, inh_syn, exc_weight, inh_weight, clump)

    """
    import itertools
    syn_list = param_list
    if not syn_list:
        syn_list = itertools.product([1, 10, 100], [1, 10, 100],  # E, I synapse numbers
                                     [1, 5, 10], [1, 5, 10],  # E, I synapse weights
                                     [1, 10, 100]  # clumping
                                     )
        syn_list = [(10, 10, 1, 1, 10), (10, 10, 5, 5, 10), (10, 10, 10, 20, 10),
                    (10, 100, 18, 6, 10),
                    (10, 10, 18, 60, 10), (10, 10*100, 18, 6, 10), (10, 10, 18, 6*100, 10)]
    load_file(filename)
    icl_rec_vec = h.Vector()
    icl_rec_vec.record(h.ldend(0.5)._ref_icl)
    h.tstop = 500
    h("forall{cli=4.25}")
    h.hoc_stdout("synapses_vs_weights.txt")
    base_gmax = 0.00005  # uS
    inh_hz = 10
    for gmax_weight in [False]:
        print("gmax_weight: {} ".format(gmax_weight))
        for (i, (exc, inh, exc_weight, inh_weight, clump)) in enumerate(syn_list):
            h.newSynapses(exc, inh)
            h.clumpSynapses(clump)
            print("clumped")
            g_rec_vec_list = list()
            for g in range(int(h.ninh)):
                g_rec_vec = h.Vector()
                g_rec_vec.record(h.GABAa[g]._ref_g)
                g_rec_vec_list.append(g_rec_vec)

            h.exHz(10000, 1, exc_weight)
            if gmax_weight:
                h.inHz(inh_hz, 1, 1)
                h.gmax_GABAa = base_gmax*inh_weight
            else:
                h.inHz(inh_hz, 1, inh_weight)
            print(syn_list[i])
            icl_mean = []
            icl_stdev = []
            spikes_mean = []
            g_total_mean = []
            for k in range(5):
                h.run()
                icl_mean.append(icl_rec_vec.mean()*1e6)
                icl_stdev.append(icl_rec_vec.stdev()*1e6)
                spikes_mean.append(h.apc.n)
                total_g = 0
                for g_rec_vec in g_rec_vec_list:
                    total_g += g_rec_vec.mean()
                g_total_mean.append(total_g)
            print("icl_mean: \tmean:\t {} \tstd:\t {}".format(np.mean(icl_mean), np.std(icl_mean)))
            print("icl_stdev: \tmean:\t {} \tstd:\t {}".format(np.mean(icl_stdev), np.std(icl_stdev)))
            print("g total: \tmean:\t {} \tstd:\t {}".format(np.mean(g_total_mean), np.std(g_total_mean)))
            print("spikes: \tmean:\t {} \tstd:\t {}".format(np.mean(spikes_mean), np.std(spikes_mean)))
