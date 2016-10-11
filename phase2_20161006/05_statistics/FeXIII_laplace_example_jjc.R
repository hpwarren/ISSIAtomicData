source("setup_v2.R")
source("laplace_functions.R")
library(parallel)  # for mclapply

npix <- 5 # try just a few pixels
nCores <- detectCores()  # This detects the # of cores available
nCores <- 4

pr_y_given_A <- matrix(NA, npix, nprior)
postmodes <- array(NA, c(npix, nprior, 2))
postcovs <- array(NA, c(npix, nprior, 4))

                                        #  This function is what is done in parallel
doPrior <- function(ii,kk){
    fn <- function(theta) log_joint_y_theta(theta, ii, Iobs[kk, ], sigmaI[kk, ],
                                            logn_grid, emissivity_grid)
    optimout <- optim(c(9.5, 9.5), fn, control = list(fnscale = -1), hessian = TRUE)
    pr_y_given_A <- laplace_approx_2D(optimout)
    postmodes <- optimout$par
    postcovs <- solve(-optimout$hessian)
    out <- list()
    out$pr_y_given_A <- pr_y_given_A
    out$postmodes <-  postmodes
    out$postcovs <-     c(postcovs)
    return(out)
}

ptm <- proc.time()

for (kk in 1:npix) {

                                        # priors done in parallel
    out <- mclapply(1:nprior, doPrior,kk = kk, mc.cores = nCores)

                                        # pull info from list
    pr_y_given_A[kk,] <- unlist(lapply(out,function(x){return(x$pr_y_given_A)}))
    postmodes[kk,, ] <- matrix(unlist(lapply(out,function(x){return(x$postmodes)})),
                               ncol = 2, byrow = TRUE)
    postcovs[kk,,] <- matrix(unlist(lapply(out,function(x){return(x$postcovs)})),
                             ncol = 4, byrow = TRUE)
    cat('Pixel ', kk, '\n')
}

print(proc.time() - ptm)

## given pixel kk and atomic data curve ii:
##   postmodes[kk, ii, 1] is the MAP estimate of log n
##   postmodes[kk, ii, 2] is the MAP estimate of log ds
##   postcovs[kk, ii, ] is var(log n), cov(log n, log ds), var(log ds)

## posterior probabilities p(A_m | y)
post_probs <- pr_y_given_A / rowSums(pr_y_given_A)

image(1:npix, 1:nprior, post_probs, ylab="Atomic Data Index (m)",
      xlab="Observation Index (k)", col=gray(seq(1, 0, length.out=64)),
      main="P(A_m | y_k)")
