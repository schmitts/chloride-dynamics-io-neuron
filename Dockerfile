#
# NEURON Dockerfile
#

# Pull base image.
FROM andrewosh/binder-base

MAINTAINER Chris Currin <chris.currin@gmail.com>

USER root

RUN \
  apt-get update && \
  apt-get install -y libncurses-dev

# update root conda environment with package requirements from environment.yml
RUN conda update -f environment.yml

# Make ~/neuron directory to hold stuff.
WORKDIR neuron

# Fetch NEURON source files, extract them, delete .tar.gz file.
RUN \
  wget http://www.neuron.yale.edu/ftp/neuron/versions/v7.7/nrn-7.7.tar.gz && \
  tar -xzf nrn-7.7.tar.gz && \
  rm nrn-7.7.tar.gz

# Fetch Interviews.
# RUN \
#  wget http://www.neuron.yale.edu/ftp/neuron/versions/v7.7/iv-19.tar.gz  && \
#  tar -xzf iv-19.tar.gz && \
#  rm iv-19.tar.gz

WORKDIR nrn-7.7

# Compile NEURON.
RUN \
  ./configure --prefix=`pwd` --without-iv --with-nrnpython=$HOME/anaconda/bin/python && \
  make && \
  make install

# Install python interface
WORKDIR src/nrnpython
RUN python setup.py install

# Add NEURON to path
ENV PATH $HOME/neuron/nrn-7.7/x86_64/bin:$PATH

# Switch back to non-root user privledges
WORKDIR $HOME
USER main
