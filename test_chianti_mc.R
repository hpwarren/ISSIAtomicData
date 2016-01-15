
#
# Instructions on loading rhdf5:
# 
#  URL: http://www.r-bloggers.com/working-with-hdf-files-in-r-example-pathfinder-sst-data/
#
# To download and compile the package:
#
# > source("http://bioconductor.org/biocLite.R")
# > biocLite("rhdf5")
#

library("rhdf5")

# read the file
fname <- "fe_13.monte_carlo.h5"
emissivity <- h5read(fname, 'emissivity')
logn <- h5read(fname, 'logn')
wavelength <- h5read(fname, 'wavelength')

# extract the dimensions
dims <- dim(emissivity)
nsim <- dims[1]
nlogn <- dims[2]
nlines <- dims[3]

# grab the first realization of the data
i <- 1
e1 <- emissivity[i,,1]
e2 <- emissivity[i,,2]
e3 <- emissivity[i,,3]
e4 <- emissivity[i,,4]
e5 <- emissivity[i,,5]

# ratio relative to the 202.044 line
ratio13 <- e1/e3
ratio23 <- e2/e3
ratio43 <- e4/e3
ratio53 <- e5/e3

# plot
par(mfrow=c(2,2))
plot(logn, ratio13, type='l')
plot(logn, ratio23, type='l')
plot(logn, ratio43, type='l')
plot(logn, ratio53, type='l')
