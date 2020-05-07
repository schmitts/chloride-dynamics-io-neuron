# coding=utf-8
"""
Set up and run a protocol.


"""
import os
import platform
import zipfile
import logging
from datetime import datetime

from neuron import h

from src.utils.file_io import cd, copy_file, create_dir
from src.utils.nrnpy import start

logging.basicConfig()
logger = logging.getLogger("run_protocol")


def run_protocol(protocol, root=None, run_using_python=True, record_output=True, timestamp=True, **kwargs):
    """
    run a test file or method
    :param protocol: the protocol to run
    :param root: name to give the run (optional)
    :param run_using_python: option for test to be processed using python or command line NEURON (default True)
    :param record_output: option to record output or not (default True)
    :param timestamp: include a timestamp for output directory
    :param kwargs: keyword arguments to be passed to test method if it is a python method
    """
    global switched_output
    # Create dir
    to_path = create_dir(root=root, timestamp=timestamp)
    date = datetime.strftime(datetime.now(), '%Y-%m-%d_%Hh%M')
    hoc_output = 'hoc_output_{}.txt'.format(date)
    print("directory: {} created".format(to_path))
    # Filter files
    included_extenstions = ['hoc', 'mod', 'py', 'm', 'txt', 'dll']
    file_names = [fn for fn in os.listdir('.')
                  if any(fn.endswith(ext) for ext in included_extenstions)]

    # Copy files
    for file_name in file_names:
        copy_file(file_name, to_path)

    print("files copied")
    # Change directory, stop if error
    with cd(to_path):
        print("changed directory\nusing environment:")
        try:
            if run_using_python:
                print('Python + NEURON')
                print("----------------------------------------------------------------------------------------")
                # load neuron standard run environment (not needed if 'from neuron import gui' used)
                h.load_file('stdrun.hoc')
                if record_output and not switched_output:
                    # to see output, open the file in an editor that can view files as they update (e.g. Sublime Text)
                    h.hoc_stdout(hoc_output)
                    switched_output = True

                if isinstance(protocol, str):
                    # protocol is a hoc file
                    print("running protocol:{}".format(protocol))
                    h.load_file(protocol)
                else:
                    print("running protocol:{}".format(protocol.__name__))
                    result = start(protocol, **kwargs)
                h.hoc_stdout()  # reset output
            else:
                print("running protocol:{}".format(protocol))
                if platform.system() == 'Windows':
                    print('Windows NEURON')
                    cmds = ["c:\\nrn\\bin\\nrniv"]
                else:
                    print('Unix NEURON')
                    cmds = ["nrniv"]
                cmds.append(protocol)
                # Record shell/terminal output (best method for large output)
                if record_output:
                    logging.basicConfig(filename=hoc_output)  # set logging to output to file
                    import tempfile
                    import subprocess
                    # no realtime output
                    # Sub_TODO: fix using http://www.cyberciti.biz/faq/python-run-external-command-and-get-output/
                    with tempfile.TemporaryFile() as tempf:
                        proc = subprocess.Popen(cmds, stdout=tempf, stderr=tempf)
                        proc.wait()
                        tempf.seek(0)
                        logger.info(tempf.read())
                else:
                    os.system(" ".join(cmds))
            print("----------------------------------------------------------------------------------------")
            print("finished running NEURON")
        except Exception as e:
            print("Exception occured: {}".format(e))
        finally:
            # Exclude specific files from being zipped and deleted:
            file_names.remove('config.hoc')
            # Zip files (not data files, log files or dlls)
            #       data (dat) and log (txt) files created AFTER file_names list was generated
            with zipfile.ZipFile('files.zip', 'w', zipfile.ZIP_DEFLATED) as zip_name:
                for file_name in file_names:
                    # exclude files from being zipped
                    if not file_name.endswith('dll'):
                        zip_name.write(file_name)
            print("code files zipped...")
            # Clean up
            for file_name in file_names:
                os.remove(file_name)
            print("and removed")

    print("done with python script")
    return to_path

