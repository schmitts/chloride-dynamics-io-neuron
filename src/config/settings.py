# coding=utf-8
"""
Global defaults, including logging format, and naming constants.
"""
import os
import platform
import logging

import colorlog
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------------------------------------------
# SET LOGGER
# ----------------------------------------------------------------------------------------------------------------------
handler = colorlog.StreamHandler()
fmt = '%(asctime)s %(name)-10s [%(filename)-10s:%(lineno)4d] %(levelname)-8s \t %(message)s'
datefmt = '%m-%d %H:%M:%S'
handler.setFormatter(colorlog.ColoredFormatter('%(log_color)s' + fmt, datefmt=datefmt))

logging.basicConfig(
        level=logging.DEBUG,
        format=fmt,
        handlers=[handler],
        datefmt=datefmt)
# to log to file, add:
# filename='/temp/myapp.log',
# filemode='w'

logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("colormath.color_conversions").setLevel(logging.WARNING)

# ----------------------------------------------------------------------------------------------------------------------
# NEURON
# ----------------------------------------------------------------------------------------------------------------------
NEURON_GUI = False
NEURON_RECOMPILE = False
NEURON_QUICK_DEFAULT = True  # use cvode where implicit
HOC_PATH = "../hoc_files/utils"
MOD_PATH = "../mod/"
NRNMECH_PATH = ''

if not os.path.isdir(MOD_PATH):
    # find the dir
    dir_path = os.path.dirname(os.path.realpath(__file__))
    MOD_PATH = os.path.join(dir_path, MOD_PATH)
    HOC_PATH = os.path.join(dir_path, HOC_PATH)
if platform.system() == 'Linux' or platform.system() == 'Darwin':
    NRNMECH_PATH = MOD_PATH + "x86_64/.libs/libnrnmech.so"
elif platform.system() == 'Windows':
    NRNMECH_PATH = MOD_PATH + "nrnmech.dll"
else:
    print("unknown system")
    exit(-1)
NRNMECH_PATH = NRNMECH_PATH.replace("\\", "/")

# ----------------------------------------------------------------------------------------------------------------------
# RANDOM
# ----------------------------------------------------------------------------------------------------------------------
RANDOM_SEED = 0
# max # events in a NetStim's stream  (Adjacent streams will be correlated by this offset.)
# // before it begins to repeat values already generated
# // by a different stream.
# // set to 0 and all NetStims will produce identical streams
RANDOM_STREAM_OFFSET = 1000

# ----------------------------------------------------------------------------------------------------------------------
# MATPLOTLIB PLOT CONFIG
# ----------------------------------------------------------------------------------------------------------------------
article_style_path = "article.mplstyle"
if not os.path.isfile(article_style_path):
    # find the file
    dir_path = os.path.dirname(os.path.realpath(__file__))
    article_style_path = os.path.join(dir_path, article_style_path)
plt.style.use(article_style_path)
logging.getLogger("settings").debug("imported style {}".format(article_style_path))

# DEFINE FIGURE USEFUL SIZES (in inches)
PAGE_W_FULL = 7.5
PAGE_H_FULL = 7.5  # make square so there's space for caption
PAGE_H_FULL_no_cap = 8.75  # no space for caption
PAGE_W_half = PAGE_W_FULL/2
PAGE_H_half = PAGE_H_FULL_no_cap/2
PAGE_W_3rd = PAGE_W_FULL/3
PAGE_H_3rd = PAGE_H_FULL_no_cap/3
PAGE_W_4th = PAGE_W_FULL/4
PAGE_H_4th = PAGE_H_FULL_no_cap/4
PAGE_W_column = 5.2  # according to https://journals.plos.org/ploscompbiol/s/figures#loc-dimensions

# DEFINE SOME HELPFUL COLOR CONSTANTS (better than normal matplotlib defaults)
default_colors = ['1f77b4', 'ff7f0e', '2ca02c', 'd62728', '9467bd', '8c564b', 'e377c2', '7f7f7f', 'bcbd22', '17becf']


class COLOR(object):
    """COLOR object for consistent choice of COLORS wherever settings.py is used
    """
    B = '#1f77b4'
    O = '#ff7f0e'
    G = '#2ca02c'
    R = '#d62728'
    Pu = '#9467bd'
    Br = '#8c564b'
    Pi = '#e377c2'
    K = '#7f7f7f'
    Ye = '#bcbd22'
    Cy = '#17becf'
    R1_B2 = '#552f72'
    R1_B3 = '#403580'
    R2_B1 = '#802456'
    R3_B1 = '#bf122b'
    E_I = '#7b4f6e'  # 50-50 mix
    # assign semantic colors
    E = R
    I = B
    A = K
    E2 = O
    NMDA = E2
    AMPA = E
    GABA = I


# ----------------------------------------------------------------------------------------------------------------------
# VARIABLE NAMES (specifically for rendering in plots)
# ----------------------------------------------------------------------------------------------------------------------

# Helper functions to clean math text
def math_clean(_s: str):
    """Remove all '$' symbols in a string"""
    return _s.replace("$", "")


def math_fix(_s: str):
    """Keep only first and last '$' in a math expression"""
    num_dollar = 0
    first = 0
    last = 0
    for idx, c in enumerate(_s):
        if c == '$':
            num_dollar += 1
            if num_dollar == 1:
                first = idx
            last = idx
    if num_dollar > 2:
        return f"{_s[:first + 1]}{math_clean(_s[first + 1:last])}{_s[last:]}"
    elif num_dollar%2 == 1:
        return f"${math_clean(_s)}$"
    else:
        return _s


# IONS
CL = cl = "$Cl\it{^-} $"
CLI = CL_i = cli = "$[{}]\mathregular{{_i}}$".format(math_clean(CL))
MILLIMOLAR = mM = "mM"
# CHLORIDE
STATIC_CHLORIDE_STR_ABBR = "Static {}".format(cl)
STATIC_CHLORIDE_STR_LONG = "Static Chloride"
DYNAMIC_CHLORIDE_STR_ABBR = "Dynamic {}".format(cl)
DYNAMIC_CHLORIDE_STR_LONG = "Dynamic Chloride"
ECL0 = f'$E{math_clean(CL)}_0$'
ECL = f'$E{math_clean(CL)}$'

TAU_KCC2 = '$\it{\\tau}_{\\rm{KCC2}}$'

# SYNAPSES
GABA = 'GABA'
EGABA = 'EGABA'
E_GABA = '$E_{GABA}$'
G_GABA = '$g_{GABA}$'
G_AMPA = '$g_{AMPA}$'
G_NMDA = '$g_{NMDA}$'
GABAA = GABAa = '$GABA_{A}$'
GABAAR = GABAaR = GABAA + 'R'
DELTA = '$\Delta $'
NABLA = '$\\nabla $'
DELTAEGABA = f'{DELTA}EGABA'
NABLAEGABA = GRADEGABA = f'{NABLA}EGABA'
