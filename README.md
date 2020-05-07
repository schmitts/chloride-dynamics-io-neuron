# chloride-dynamics-io-neuron
Code for the paper Chloride dynamics alter the input-output properties of neurons

Authors: 
- [Christopher Currin](https://chriscurrin.com)

## 1. Create environment from requirements using [conda](https://docs.conda.io/en/latest/)

`conda env create -f environment.yaml`

Target Python version: >= 3.6

> development started on Python 2.7, so there may be *some* legacy code. 

## 2. Install [NEURON + Python](https://www.neuron.yale.edu/neuron/)  

## 3. Play around with notebook

`jupter notebook opdynamics.ipynb`

if the python kernel for the created environment (from step 1) is not showing up:

```bash 
conda activate chloride-dynamics-io-neuron 
python -m ipykernel install --user --name chloride-dynamics-io-neuron
```

# Notes
- NMODL (`.mod`) files are automatically attempted to be compiled in `shared.INIT_NEURON()` on any system.
  \
  For those on Windows, I've included a `.bat` file in the `mod` folder to help compile the files if you want to run things independently. 
- A lot of initial work was done in NEURON's hoc language, so Python here is used as a scripting language mostly.
- sorry about the MATLAB file, I was young...
