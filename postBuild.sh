#!/usr/bin/env bash

# Fetch NEURON source files, extract them, delete .tar.gz file.
wget http://www.neuron.yale.edu/ftp/neuron/versions/v7.7/nrn-7.7.tar.gz && \
  tar -xzf nrn-7.7.tar.gz && \
  rm nrn-7.7.tar.gz

# Fetch Interviews.
# wget http://www.neuron.yale.edu/ftp/neuron/versions/v7.7/iv-19.tar.gz  && \
#  tar -xzf iv-19.tar.gz && \
#  rm iv-19.tar.gz

echo "compiling NEURON"

cd nrn-7.7

./configure --prefix=`pwd` --without-iv --with-nrnpython=/srv/conda/envs/notebook/bin/python && \
  make && \
  make install

echo "installing NEURON+Python"
cd src/nrnpython

python setup.py install

PATH=$HOME/nrn-7.7/x86_64/bin:$PATH

