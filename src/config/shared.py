# coding=utf-8
from __future__ import print_function, division

import glob
import hashlib
import logging
import numpy as np
from neuron import h

from src.config import settings

logger = logging.getLogger('shared')

t_vec = None
initialised = False
__KWARGS__ = {}


def INIT_NEURON(reinit=False):
    """ Initialise NEURON

    - Compile mod files in `settings.MOD_PATH`. Files are recompiled if there are changes or
        `settings.NEURON_RECOMPILE` is True.
    - Load mod files into neuron


    :param reinit:
    :type reinit:
    :return:
    :rtype:
    """
    global initialised, t_vec, __KWARGS__
    if initialised and not reinit:
        return True
    else:
        initialised = True
    # compile mod files
    if __mod_files_changed(settings.MOD_PATH) or settings.NEURON_RECOMPILE:
        from src.utils.compile_mod import compile_mod
        output = compile_mod(path=settings.MOD_PATH, mod=True)
        if "Error" in str(output):
            raise Exception("MOD FILES not compiled successfully")
    # load mod files
    h.nrn_load_dll(settings.NRNMECH_PATH)
    # load hoc files including usefulFns.hoc
    for hoc_file in glob.glob(settings.HOC_PATH + "/*.hoc"):
        h.load_file(hoc_file.replace("\\", "/"))
    # show GUI
    if settings.NEURON_GUI:
        # noinspection PyUnresolvedReferences
        from neuron import gui
        h.showRunControl()

    # general properties
    h.celsius = 37
    h.v_init = -65
    h.random_stream_offset_ = settings.RANDOM_STREAM_OFFSET
    logger.info("celsius={} and v_init={}".format(h.celsius, h.v_init))
    t_vec = h.Vector()
    np.random.seed(settings.RANDOM_SEED)

    # reset __KWARGS__
    __KWARGS__ = {}
    env_var(celsius=h.celsius, v_init=h.v_init)


def env_var(**kwargs):
    """Share environment variables (the NEURON kind, not the OS kind), by using this method to assign and call to
    retrieve."""
    if kwargs:
        for k, v in kwargs.items():
            __KWARGS__[k] = v
    return __KWARGS__


def show_n(n=1, ion=False):
    """Show n figures, closes the rest."""
    import matplotlib.pyplot as plt
    for i in plt.get_fignums():
        if i > n:
            plt.close(i)
    if ion:
        plt.ion()
    plt.show()


def hashable(cls):
    """Decorator to make a hashable class by using its `__str__` method.

    The class can thus be customised to be hashed and cached based on what properties are in its `__str__`
    """
    def __hash__(self):
        return hashlib.md5(str(self).encode('utf-8')).hexdigest()

    def hash_extra(self, extra=""):
        full_str = str(self) + extra
        return hashlib.md5(full_str.encode('utf-8')).hexdigest()

    cls.__hash__ = __hash__
    cls.hash_extra = hash_extra
    return cls


def __mod_files_changed(path=settings.MOD_PATH):
    md5_files = glob.glob(path + "/hash.md5")
    if len(md5_files) == 0:
        old_md5 = ''
    elif len(md5_files) == 1:
        with open(md5_files[0]) as f:
            old_md5 = f.read()
    else:
        raise BaseException("Too many hash files")

    new_md5_hash = hashlib.md5()
    for mod_file in glob.glob(path + "/*.mod"):
        new_md5_hash.update(__md5(mod_file).encode('utf-8'))
    new_md5 = new_md5_hash.hexdigest()
    if new_md5 == old_md5:
        return False
    else:
        # there are changes
        with open(path + "/hash.md5", 'w') as hash_file:
            hash_file.write(new_md5)
        logger.info("there were changes in the mod file directory")
        return True


def __md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(2 ** 20), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()
