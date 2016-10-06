require(rhdf5)

# read the file
fname <- "../02_test/test_intensities_fe_13.h5"

emissivity_grid <- h5read(fname, "emissivity")
logn_grid <- h5read(fname, "logn")
wavelength <- h5read(fname, "wavelength")
Iobs <- h5read(fname, "intensities")
sigmaI <- h5read(fname, "intensities_error")

dims <- dim(emissivity_grid)
nprior <- dims[1]
ngrid <- dims[2]
nlines <- dims[3]

print(nprior)
print(ngrid)
print(nlines)
