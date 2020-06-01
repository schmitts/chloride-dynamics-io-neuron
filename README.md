# [Chloride dynamics alter the input-output properties of neurons](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007932)

Currin CB, Trevelyan AJ, Akerman CJ, Raimondo JV (2020) Chloride dynamics alter the input-output properties of neurons. PLOS Computational Biology 16(5): e1007932. https://doi.org/10.1371/journal.pcbi.1007932

Code authors: 
- [Christopher Currin](https://chriscurrin.com)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ChrisCurrin/chloride-dynamics-io-neuron/master?filepath=chloride-dynamics-io-neuron.ipynb)


# To run on your own computer
## 1. Create environment from requirements using [conda](https://docs.conda.io/en/latest/)

`conda env create -f environment.yaml`

Target Python version: >= 3.6

> development started on Python 2.7, so there may be *some* legacy code. 

> `environment.yml` includes installation for [NEURON + Python](https://www.neuron.yale.edu/neuron/) via the `neuron` package.  
> feel free to remove this from `environment.yml` if it may conflict with your local version

## 2. Play around with notebook

`jupter notebook chloride-dynamics-io-neuron.ipynb`

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
