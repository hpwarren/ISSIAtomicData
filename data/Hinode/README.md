
#### Background

Let's assume that the intensity of each line can be described by
\begin{equation}
  I_\lambda = \epsilon_\lambda(n_e, T_e)\,n_e^2\,ds
\end{equation}

where $\epsilon_\lambda$ is the atomic data for that line, $n_e$ is the electron density, $T_e$ is
the electron temperature, and $ds$ is the path length through the solar atmosphere. We assume that
all of the Fe XIII emission is formed at the same temperature and use a fixed $T_e$.

#### Analyzing the Observations

One approach to utilizing the observed intensities is to use equation 1 to form ratios and to infer
densities from the observed and theoretical ratios. This is what we've discussed previously.

Another approach is to find the electron density and path length that fits the observed intensities
as closely as possible. I've written a routine that implements this. Here is an example fit using
the default CHIANTI atomic data:

```
  model log_n = 9.48 +- 0.016
 model log_ds = 9.14 +- 0.031
         chi2 = 153.2
      Line    Imodel      Iobs    SigmaI  dI/Sigma      dI/I
   200.021   2064.02   1809.67     32.90      7.73      14.1
   201.121   2506.85   2946.72     51.14      8.60      14.9
   202.044   4299.64   4153.84     64.85      2.25       3.5
   203.152    932.87   1071.74     48.20      2.88      13.0
   203.826  10223.57  10620.57    160.95      2.47       3.7
```

Again, the idea is to fit all of the Fe XIII intensities as closely as possible simultaneously. For
these parameters the spectrum looks something like this

![A spectrum computed from CHIANTI for the best-fit density and path length derived from the Fe
 XIII intensities.](fe_13_fit_intensities_spec.jpg)

We can repeat this calculation using the perturbed CHIANTI atomic data. For the perturbed data we
also appear to get reasonable fits to the observed spectra. For example, 

```
  model log_n = 9.80 +- 0.026
 model log_ds = 8.53 +- 0.050
         chi2 = 68.5
      Line    Imodel      Iobs    SigmaI  dI/Sigma      dI/I
   200.021   1775.97   1809.67     32.90      1.02       1.9
   201.121   2608.44   2946.72     51.14      6.62      11.5
   202.044   4283.36   4153.84     64.85      2.00       3.1
   203.152    997.34   1071.74     48.20      1.54       6.9
   203.826  11290.98  10620.57    160.95      4.17       6.3
```

We can repeat this for all of the available perturbed atomic data, which have been calculated for
5-30%. In the following figure we see two features. On the left we see that the path
length and the density are inversely correlated, as we expect from Equation 1. On the right we see
that for many realizations of the atomic data there are parameter sets that fit the observations as
well or better than the default parameters.

![A summary of 600 fits to the observations using perturbed atomic data. The important part of this
 figure is the upper right, which shows the best-fit density as a function of chi-squared. It
 suggests that strongly perturbed version of CHIANTI could also reproduce the observed
 intensities. The bottom panels show the distributions of density and path length for chi-squared
 less than 250. The large dot represents the solution with the default version of
 CHIANTI.](fe_13_fit_intensities_plot.jpg)