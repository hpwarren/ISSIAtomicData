interpolate_emissivities <- function(logn, index, logn_grid, emissivity_grid) {
  emis <- rep(NA, 5)
  for (ii in 1:nlines) {
    emis[ii] <- approx(logn_grid, emissivity_grid[index, , ii], logn)$y
  }
  emis
}

log_joint_y_theta <- function(theta, index, Iobs, sigmaI, logn_grid, 
                              emissivity_grid)
{
  logn <- theta[1]
  logds <- theta[2]

  A <- interpolate_emissivities(logn, index, logn_grid, emissivity_grid)

  dunif(logn, 7, 12, log=TRUE) + dcauchy(logds, 9, 5, log=TRUE) + # log prior
    sum(dnorm(Iobs, A * 10^(2 * logn + logds), sigmaI, log=TRUE)) # log lik
}

laplace_approx_2D <- function(optimout) {
  2 * pi / sqrt(det(-optimout$hessian)) * exp(optimout$value)
}
