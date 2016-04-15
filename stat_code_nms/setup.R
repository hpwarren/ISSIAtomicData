require(rhdf5)

# read the file
fname <- "../atomic_data_v2_gdz/fe_13.monte_carlo_normal_nsim=1000.h5"
emissivity_grid <- h5read(fname, "emissivity")
logn_grid <- h5read(fname, "logn")
wavelength <- h5read(fname, "wavelength")

dims <- dim(emissivity_grid)
nprior <- dims[1]
ngrid <- dims[2]
nlines <- dims[3]

# read the scaling factor
# this term multiplies all of the emissivities
fname_cf <- "../data/HinodeAnalysisGDZ_v1/chianti_factor.fe13.h5"
cf <- h5read(fname_cf, "chianti_factor")

for (ii in 1:nprior) {
  for (jj in 1:nlines) {
    emissivity_grid[ii, , jj] = emissivity_grid[ii, , jj] * cf
  }
}

# read the observed intensities
fname_eis <- "../data/HinodeAnalysisGDZ_v1/eis_l1_20130708_002042.fe_density.h5"
eis_ints <- h5read(fname_eis, 'intensities')
eis_err <- h5read(fname_eis, 'intensities_error')
eis_wave <- h5read(fname_eis, 'wavelength')


# use index to access the observed eis intensities
index <- c(1,2,4,5,3)

Iobs <- eis_ints[, index]
sigmaI <- eis_err[, index]

