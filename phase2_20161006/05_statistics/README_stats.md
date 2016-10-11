#### ISSIAtomicData/phase2_20161006/05_statistics

The command that you want:

```
 source("FeXIII_laplace_example.R")
```

or

```
 source("FeXIII_laplace_example_jjc.R")
```

for Jessi's parallel version.


These R routines illustrate how to compute the models for the observed and test intensities.  The
setup file looks like this:

```
require(rhdf5)

## read the intensity file of interest
f = c("../04_observed/eis_l1_20130708_002042.fe_density.h5",
      "../02_test/test_intensities_fe_13.h5")
print(f)
res <- readline("which file do you want to process? ")
res <- as.numeric(res)
fname <- f[res]
print(fname)

emissivity_grid <- h5read(fname, "emissivity")
logn_grid <- h5read(fname, "logn")
wavelength <- h5read(fname, "wavelength")
Iobs <- h5read(fname, "intensities")
sigmaI <- h5read(fname, "intensities_error")

dims <- dim(emissivity_grid)
nprior <- dims[1]
ngrid <- dims[2]
nlines <- dims[3]

if (fname == "../02_test/test_intensities_fe_13.h5") {
    ## these variables are defined for the test data, where logn and logds are known
    logn_model <- h5read(fname, "logn_obs")
    logds_model <- log10( h5read(fname, "ds_obs") )
} else {
    ## these variables are undefined for the real data
    logds_model = NULL
    logn_model = NULL
}
```