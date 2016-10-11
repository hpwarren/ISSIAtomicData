source("setup_v2.R")
source("laplace_functions.R")

npix <- 25 # try just a few pixels

pr_y_given_A <- matrix(NA, npix, nprior)
postmodes <- array(NA, c(npix, nprior, 2))
postcovs <- array(NA, c(npix, nprior, 2, 2))

ptm <- proc.time()

for (kk in 1:npix) {
    for (ii in 1:nprior) {
        fn <- function(theta) log_joint_y_theta(theta, ii, Iobs[kk, ], sigmaI[kk, ],
                                                logn_grid, emissivity_grid)
        optimout <- optim(c(9, 9), fn, control = list(fnscale = -1), hessian = TRUE)
        pr_y_given_A[kk, ii] <- laplace_approx_2D(optimout)
        postmodes[kk, ii, ] <- optimout$par
        postcovs[kk, ii, , ] <- solve(-optimout$hessian)

        if (ii == 1 || ii %% 25 == 0){
            if(!is.null(logds_model) && !is.null(logn_model)) {
                this_logds_model = logds_model[kk]
                this_logn_model = logn_model[kk]
            } else {
                this_logds_model = NULL
                this_logn_model = NULL
            }
            print_results(postmodes[kk,ii,1], postmodes[kk,ii,2], ii, logn_grid, emissivity_grid,
                          Iobs[kk,], sigmaI[kk,], logn_model=this_logn_model, 
                          logds_model=this_logds_model)
        }
    }
    mean_logn = mean(postmodes[kk,,1])
    mean_logds = mean(postmodes[kk,,2])
}

print(proc.time() - ptm)

## given pixel kk and atomic data curve ii:
##   postmodes[kk, ii, 1] is the MAP estimate of log n
##   postmodes[kk, ii, 2] is the MAP estimate of log ds
##   postcovs[kk, ii, , ] is the estimated covariance matrix of (log n, log ds)

## posterior probabilities p(A_m | y)
post_probs <- pr_y_given_A / rowSums(pr_y_given_A)

image(1:npix, 1:nprior, post_probs, ylab="Atomic Data Index (m)",
      xlab="Observation Index (k)", col=gray(seq(1, 0, length.out=64)),
      main="P(A_m | y_k)")
