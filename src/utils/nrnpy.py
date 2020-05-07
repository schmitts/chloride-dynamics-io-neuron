# coding=utf-8
"""
NEURON-related python methods
"""
import os

from neuron import h
import logging

logger = logging.getLogger("nrnpy")


def load_file(filename, set_file_name=True, force_reload=True):
    import glob
    if filename.endswith(".hoc"):
        filename = filename[0:-4]
    src = 'src/' if os.path.exists('src') else ''
    file_path = glob.glob(f"{src}hoc_files/**/{filename}.hoc", recursive=True)[0].replace('\\', '/')
    h.load_file(int(force_reload), file_path)
    if set_file_name:
        h("strdef fileName")
        h.fileName = filename
    logger.info("opened: {}.hoc \t reloaded={}".format(filename, force_reload))


def start(method, filenames=None, param_list=None, config_functions=None):
    """
    start a Python+NEURON method as a protocol
    :param method: python method to run
    :param filenames: list of files to protocol
    :param param_list: list of parameters for method (optional). See defaults in passed method
    :param config_functions: the configurations to (over) load, and the args - {function_name:{"arg1":"arg1value"}}
    :return: 1 if no filenames specified, 0 if all filenames successfully tested using method
    """
    if filenames is None:
        print("no filenames specified")
        return 1
    h.load_file("usefulFns.hoc")
    h.load_file("cell.hoc")
    if config_functions is not None:
        for config_function, kwargs in config_functions.items():
            if kwargs is None:
                config_function()
            else:
                config_function(**kwargs)
    h("strdef fileName, saveName, openFile")
    for filename in filenames:
        method(filename=filename, param_list=param_list)
    else:
        return 0


def save_run(savename, set_cli=None, skip_if_exists=False):
    """

    :param savename:
    :param set_cli:
    """
    h.saveName = savename
    if skip_if_exists:
        from os import listdir
        from os.path import isfile, join
        onlyfiles = [f for f in listdir(path='..') if isfile(join('.', f))]
        for f in onlyfiles:
            if f.startswith(savename):
                print("skipping {} because it has already been run".format(h.saveName))
                return
    print("running {} and saving {} ...".format(h.fileName, h.saveName))
    # load custom run file
    h.load_file("run.hoc")
    # doRun() is procedure in run.hoc
    if set_cli is not None:
        h.doRun(set_cli)
    else:
        h.doRun()


def get_base_vm_cli(file_name, load=False, compartment=None):
    """
    Determine steady-state internal chloride concentration by
    1) instantiating a class that extends BaseNeuron (NOTE: class should not have spiking at default values)
    2) adding KCC2 using the add_kcc2() method
    3) setting [Cl]_i to an arbitrary value (can be set in method invocation)
    running a fast simulation for a long time
    checking if chloride is at steady state (repeat 4 until it is)
    return steady state Vm and [Cl]_i

    :return: steady state Vm and [Cl]_i
    """
    logger.debug(f"finding steady-state vm and cli for {file_name}")
    if load:
        load_file(file_name, set_file_name=False)
    if compartment is None:
        compartment = h.soma
    h.inPy(0)
    h.ex(0)
    if file_name.find("KCC2") > 0:
        kcc2_already = True
    else:
        kcc2_already = False
        h.load_file(1, "insert_KCC2.hoc")
        logger.info("KCC2 temporarily inserted")

    h.tstop = 50000
    h.useCV()
    h.finitialize(-65)

    h.run()

    def at_steady_state(compartment, continue_dt):
        """
        check if [Cl]_i is at steady state
        :param continue_dt: amount of time to run
        :return: [Cl]_i if at steady state, False otherwise
        """
        v_start = compartment(.5).v
        cli_start = compartment(.5).cli
        h.continuerun(h.tstop + continue_dt)
        h.tstop += continue_dt
        v_after = compartment(.5).v
        cli_after = compartment(.5).cli
        if v_after - v_start < 1e-8 and cli_after - cli_start < 1e-8:
            return cli_after
        else:
            return False

    num_steady_state_checks = 0
    while not at_steady_state(compartment, 1):
        h.continuerun(h.tstop + 10000)
        h.tstop += 10000
        num_steady_state_checks += 1
        if num_steady_state_checks > 10:
            logger.info("not at steady state even after {} ms".format(50000 + num_steady_state_checks*10000))
            exit(-1)

    h.disableCV()
    vm, cli = compartment(.5).v, compartment(.5).cli
    logger.info("steady state [Cl]_i {}".format(cli))
    logger.info("steady state Vm {}".format(vm))
    logger.info("took {} ms (simulation time)".format(50000 + num_steady_state_checks*10000))
    if not kcc2_already:
        h("""
            forall{
                uninsert KCC2
            }
        """)
        logger.info("temporary KCC2 removed")
    return vm, cli


def save(title, t_vec, v_vec):
    """Save NEURON voltage recoding to file (with time points)."""
    import pickle
    with open(title + '_t_vec.p', 'w') as t_vec_file:
        pickle.dump(t_vec.to_python(), t_vec_file)
    with open(title + '_v_vec.p', 'w') as v_vec_file:
        pickle.dump(v_vec.to_python(), v_vec_file)


def load(title):
    """Load NEURON voltage recording from file (with time points)."""
    import pickle
    # Unpickle
    with open(title + '_t_vec.p') as t_vec_file:
        py_t_vec = pickle.load(t_vec_file)
    t_vec_restore = h.Vector(py_t_vec)
    with open(title + '_v_vec.p') as vec_file:
        py_v_vec = pickle.load(vec_file)
    v_vec_restore = h.Vector(py_v_vec)
    return t_vec_restore, v_vec_restore
