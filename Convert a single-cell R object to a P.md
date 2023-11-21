## Convert a single-cell R object to a Python anndata object or vice-versa

https://github.com/JiekaiLab/scDIOR

# scDior converts an R SCE object to a .h5 file.
# Then, use the helper functions from diopy (the python version of scDIOR) to import the first time. Then, you can save it as a .h5ad

#
adata = diopy.input.read_h5(file = 'scdata.h5')
