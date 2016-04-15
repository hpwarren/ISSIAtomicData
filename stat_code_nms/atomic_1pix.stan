functions {
  real interpolate(real x_new, vector x, vector y) {
    int ii;
    real run;
    real wt;
    int Nx;

    Nx <- num_elements(x);
    ii <- 2;
    run <- x[ii];

    // inefficient but straightforward...
    while (x_new > run) {
      ii <- ii + 1;
      run <- x[ii];
    }
    wt <- (x_new - x[ii - 1]) / (x[ii] - x[ii - 1]);
    return (1 - wt) * y[ii - 1] + wt * y[ii];
  }
}
data {
  int<lower=1> nlines;
  int<lower=1> nprior; // number of samples from emissivity prior distribution
  int<lower=1> ngrid; // number of grid points for emissivity
  vector[ngrid] logn_grid;
  real<lower=0> emissivity_grid[nprior, ngrid, nlines];
  vector[nlines] Iobs; // observed intensity
  vector[nlines] sigmaI; // standard deviation on observed intensity
}
transformed data {
  real minlogn;
  real maxlogn;
  vector[ngrid] emissivity_vec_arr[nprior, nlines];

  minlogn <- min(logn_grid);
  maxlogn <- max(logn_grid);

  // rearrange emissivity grid as array of vectors
  // for faster indexing when interpolating
  for (ii in 1:nprior) {
    for (jj in 1:nlines) {
      for (kk in 1:ngrid) {
        emissivity_vec_arr[ii, jj, kk] <- emissivity_grid[ii, kk, jj];
      }
    }
  }
}
parameters {
  real<lower=minlogn, upper=maxlogn> logn;
  real logds;
}
transformed parameters {
  vector[nprior] log_summand_lik; // log summand of p(Iobs | logn, logds)
  
  log_summand_lik <- rep_vector(0, nprior);
  for (ii in 1:nprior) {
    for (jj in 1:nlines) {
      log_summand_lik[ii] <- log_summand_lik[ii] + 
        normal_log(Iobs[jj], 
          interpolate(logn, logn_grid, emissivity_vec_arr[ii, jj]) * 
            pow(10, 2 * logn + logds),
          sigmaI[jj]);
    }
  }
}
model {
  // implicitly: logn ~ uniform(minlogn, maxlogn);
  logds ~ cauchy(9, 5);
  increment_log_prob(log_sum_exp(log_summand_lik));
}
generated quantities {
  // prob of emissivity curve given (Iobs, logn, logds)
  vector<lower=0, upper=1>[nprior] condprob;
  condprob <- softmax(log_summand_lik);
}
