The latest version can be installed using pip

## Installing via pip

`pip install transit-chem`

In some instances (particularly windows machines), you may run into
trouble while trying to install [numba](http://numba.pydata.org/). 
If this happens, The most reliable way I've found is to use 
[conda](https://docs.conda.io/en/latest/) to install numba. 


## Installing in a conda environment

```shell
conda create -n transit-chem python=3.7
conda activate transit-chem
conda install numba
pip install transit-chem
```

