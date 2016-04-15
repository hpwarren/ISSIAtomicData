# read in observations and atomic data curves
source("setup.R")

# Stan is available at: mc-stan.org
# Installation instructions for Stan in R: 
# https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
library(rstan)

pixel_index <- 661

atomic_dat <- list(nlines = nlines,
                   nprior = nprior,
                   ngrid = ngrid,
                   logn_grid = logn_grid,
                   emissivity_grid = emissivity_grid,
                   Iobs = Iobs[pixel_index, ],
                   sigmaI = sigmaI[pixel_index, ])

# This will take some time to run:
atomic_stanfit <- stan(file = "atomic_1pix.stan", 
                       data = atomic_dat, 
                       iter = 4000, 
                       chains = 5)

print(atomic_stanfit, pars = c("logn", "logds"))

traceplot(atomic_stanfit, pars = c("logn", "logds"))

atomic_fit_extracted <- extract(atomic_stanfit)

# posterior probabilities p(A_m | y)
post_probs <- colMeans(atomic_fit_extracted$condprob)

