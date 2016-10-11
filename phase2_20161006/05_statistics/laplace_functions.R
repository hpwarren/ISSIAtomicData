interpolate_emissivities <- function(logn, index, logn_grid, emissivity_grid) {
    emis <- rep(NA, nlines)
    for (ii in 1:nlines) {
        emis[ii] <- approx(logn_grid, emissivity_grid[index, , ii], logn)$y
    }
    emis
}

print_results <- function(logn, logds, index, logn_grid, emissivity_grid, Iobs, sigmaI,
                          logn_model=NULL, logds_model=NULL){

    A <- interpolate_emissivities(logn[1], index, logn_grid, emissivity_grid)
    ints_model = A * 10^(2 * logn + logds)

    if(!is.null(logds_model) && !is.null(logn_model)) {
        message( sprintf("  logn = %6.2f [input %5.2f] ", logn, logn_model) )
        message( sprintf(" logds = %6.2f [input %5 .2f] ", logds, logds_model) )
        
    } else {
        message( sprintf("  logn = %10.2f", logn) )
        message( sprintf(" logds = %10.2f", logds) )
    }
    message( sprintf(" %8s %8s  %8s ", "Obs", "Model", "Ratio") )
    for (ii in 1:nlines){
        message( sprintf(" %8.2f  %8.2f   r = %4.2f", Iobs[ii], ints_model[ii],
                         ints_model[ii]/Iobs[ii]) )
    }

}

log_joint_y_theta <- function(theta, index, Iobs, sigmaI, logn_grid, emissivity_grid)
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
