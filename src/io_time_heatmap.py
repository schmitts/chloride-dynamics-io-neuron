# coding=utf-8
from src.config import shared

import os
import time
from matplotlib import pyplot as plt

from src.iocurves.run_range import run_precise
from src.iocurves.sim import show_menu

from src.iocurves.vis import plot_io_curve, plot_compare_dcl

shared.INIT_NEURON()

# fix python only using one core
# (https://stackoverflow.com/questions/15639779/why-does-multiprocessing-use-only-a-single-core-after-i-import-numpy)
os.system("taskset -p 0xff %d"%os.getpid())
fig_time = time.strftime('%Y_%m_%d_%Hh%Mm')
is_menu = 0
save_location = "results_plots/"


# noinspection PyArgumentEqualDefault
def iotime(file_name="distal",
           ifr_windowsize=0.02,
           time_points=None,
           synapse_type=1,  # 1 for persistent
           synapse_numbers=(100, 100),
           inh_input=None,
           exc_input=None,
           diam=None,
           trials=5,
           save=True):
    """Calculate and plot heatmap of input-output curves for persistent synapses on a distal dendrite."""
    if exc_input is None:
        exc_input = list(range(0, 33, 3))
    if inh_input is None:
        inh_input = list(range(0, 22, 2))
    if time_points is None:
        time_points = [int(ifr_windowsize*1000), 100, 500, 1000, 10000]
    if diam is None:
        diam = {'ldend': 1}
    show_menu()
    ifr_windowsize = 0.02  # in seconds
    time_points = [int(ifr_windowsize*1000), 100, 500, 1000, 10000]
    save_formats = ['png', 'eps']
    common_args = dict(file_name=file_name,
                       plot=False,
                       ifr_windowsize=ifr_windowsize,
                       synapse_type=synapse_type,
                       synapse_numbers=synapse_numbers,
                       diam=diam,
                       trials=trials,
                       save=save)
    # static_runs, kcc2_runs = run_dynamic(0, 20, 0, 20, mid='half', fr_diff=10, **common_args)
    static_runs, kcc2_runs = run_precise(inh_input=inh_input, exc_input=exc_input, **common_args)
    static_runs_keys = static_runs.keys()
    static_runs_keys = []
    for static_name in static_runs_keys:
        file_name = "{save_location}traces_{result_name}.{format}".format(save_location=save_location,
                                                                          result_name=static_name, format="eps")
        if os.path.exists(file_name):
            continue
        static_run = static_runs[static_name]
        underscore = static_name.find("_")
        kcc2_run = kcc2_runs[static_name[:underscore] + "_KCC2" + static_name[underscore:]]
        f, axes = plot_compare_dcl(static_name, [static_run, kcc2_run],
                                   ifr_windowsize=ifr_windowsize, combine=False)
        # f.savefig(file_name)

    f, axes = plot_io_curve([static_runs, kcc2_runs], ifr_windowsize=ifr_windowsize, combine=True, heatmap=True,
                            fill=True,
                            time_points=time_points, show_cli=True,
                            save_args=dict(fig_name="dynamic", formats=save_formats))
    # f, axes = plot_io_curve([static_runs, kcc2_runs], ifr_windowsize=ifr_windowsize, combine=True, heatmap=True,
    #                    time_points=None,
    #                    save_args=dict(fig_name="dynamic", formats=save_formats))
    return static_runs, kcc2_runs


if __name__ == '__main__':
    iotime()
    plt.show()
