
library("rhdf5")

# a function to compute the model intensities using the CHIANTI emissivities

fe_13_compute_intensities <- function(logn, emissivity, input_logn, input_logds) {

  n_lines <- dim(emissivity)[2]
  ints_model <- array(n_lines)

  for (i in 1:n_lines)
  {
    res <- approx(logn, emissivity[,i], input_logn)
    this_emiss = res[[2]]
    em <- 10.0^(2*input_logn + input_logds)
    ints_model[i] <- this_emiss*em
  }

  return(ints_model)
}

# read the perturbed atomic data file, first element is the default
fname_rnd <- "../../atomic_data_v2_gdz/fe_13.monte_carlo_normal_nsim=1000.h5"
emissivity <- h5read(fname_rnd, 'emissivity')
logn <- h5read(fname_rnd, 'logn')
wavelength <- h5read(fname_rnd, 'wavelength')
n_lines = dim(emissivity)[3]
n_emiss = dim(emissivity)[1]

# read the scaling factor that should have been included in the beginning!
# this term multiplies all of the emissivities
fname_cf <- "../HinodeAnalysisGDZ_v1/chianti_factor.fe13.h5"
cf <- h5read(fname_cf, 'chianti_factor')

for (i in 1:n_emiss)
{
	for (j in 1:n_lines)
	{
		emissivity[i,,j] = emissivity[i,,j]*cf
  }
}

# read the observed intensities
fname_eis <- "../HinodeAnalysisGDZ_v1/eis_l1_20130708_002042.fe_density.h5"
eis_ints <- h5read(fname_eis, 'intensities')
eis_err <- h5read(fname_eis, 'intensities_error')
eis_wave <- h5read(fname_eis, 'wavelength')

# another small problem, the observed intensities and the emsisivites are not ordered the same way!
# use index to access the observed eis intensities
index <- c(1,2,4,5,3)
print("Wavelengh Check")
print(eis_wave[index])
print(wavelength)

p <- 661
this_eis_ints <- eis_ints[p,index]
this_eis_err <- eis_err[p,index]

# read the fit parameters that we've already calculated for this pixel
fname_param = "../HinodeAnalysisGDZ_v1/fe_13_fit_intensities.normal.h5"
param = h5read(fname_param, 'res')

# --- show the default result, #1

this_param = param[1,]
ints_model <- fe_13_compute_intensities(logn, emissivity[1,,], this_param[1], this_param[2])

print(sprintf(" log_n = %6.2f", this_param[1]))
print(sprintf("log_ds = %6.2f", this_param[2]))
print(sprintf("%10s %10s %10s %10s", "WAVE", "MODEL", "OBS", "ERR"))
for (j in 1:n_lines)
{
  s = sprintf("%10.3f %10.2f %10.2f %10.2f", wavelength[j], ints_model[j], this_eis_ints[j],
		this_eis_err[j])
  print(s)
}
