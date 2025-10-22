Code for Shuwen Tan's NLIW model analysis.


For development, install the package in a new conda environment as follows:
```
conda env create -f environment.yml
conda activate NLIW
```

Then install xarray
```
conda install -c conda-forge xarray dask netCDF4 bottleneck
```

to self: need to use pip within environment to install any other packages
/Users/stan/miniforge3/envs/NLIW/bin/pip install -e .
/Users/stan/miniforge3/envs/NLIW/bin/pip install xarray
/Users/stan/miniforge3/envs/NLIW/bin/pip install -U scikit-image

check version, update pip, and install scikit-image:
```
python -c "import skimage; print(skimage.__version__)"
python -m pip install -U pip
python -m pip install -U scikit-image
python -c "import xarray; print(xarray.__version__)"
python -m pip install -U xarray
pip install pyproj
pip install cmocean
```


If developing with jupyter notebook/lab it may be convenient to install the kernel as follows:
```
conda activate NLIW
python -m ipykernel install --user --name NLIW --display-name NLIW
```