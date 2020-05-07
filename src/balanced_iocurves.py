# coding=utf-8
"""
Run and plot balanced input-output curves.

Note that figure 3 in the paper was made in Matlab.
"""
import glob
import os
import logging
import pandas as pd
import matplotlib.pyplot as plt

from src.iocurves.balanced_input.vis import plot_balanced_iocurves
from src.protocols import protocol_balanced_synaptic_input
from src.utils.file_io import read_hoc_output

logging.basicConfig(level=logging.DEBUG)

logger = logging.getLogger("balanced input")


def balanced_io_curves():
    """
    Run and plot balanced input-output curves.
    """
    prox = [(260, 30),
            (290, 60),
            (330, 90),
            (490, 180),
            (780, 330),
            ]
    distal = [(200, 30),
              (230, 140),
              (250, 300),
              (270, 500),
              (300, 800),
              ]
    data_dirs = protocol_balanced_synaptic_input(distal, prox)
    for directory in data_dirs:
        data_folder = glob.glob(os.path.join(directory, "hoc_output*"))
        data_path = os.path.join(directory, "df.h5")
        if os.path.exists(data_path):
            df = pd.read_hdf(data_path, key='df')
        else:
            df = read_hoc_output(data_folder)
            df.to_hdf(data_path, key="df")

        if 'proximal' in data_path:
            cmap = 'Greens'
            ylims = (0, 90)
        else:
            cmap = 'Blues'
            ylims = (0, 120)

        ax = plot_balanced_iocurves(df, cmap)
        ax.set_ylim(*ylims)

    return ax


if __name__ == '__main__':
    balanced_io_curves()
    plt.show()
