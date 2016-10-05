
#### ISSIAtomicData/phase2_20161006/02_test

* TEST_INTENSITIES_FE_13: This routine generates fake intensities using an assumed electron density
  and path length. Some nuances

    - The emissivities computed 01_chianti_errors need some 'extra' factors (e.g., elemental
      abundance, ionization fraction) to be useful for comparison with observations in absolute
      units. They are added in here.

    - We select random densities uniformly on the range 8.5 to 11.0.

    - We estimate the path length from an approximation derived from a steady, uniform heating model
  	  (Martens et al. 2000, equation 24).
      \begin{equation}
        ds = \frac{2.56\times10^8}{P_0}
		  \end{equation}
      where $P_0=2k_bn_eT_e$ is in dyne cm$^{-2}$ and $ds$ is in cm. The intensity is
			\begin{equation}
			 I_\lambda = \epsilon_\lambda(n_e, T_e) n_e^2\,ds
		 \end{equation}
		 where we use the default CHIANTI emissivity to compute the intensity.

    - Statistical uncertatinties are added assuming the EIS pre-flight effective areas, a $60\,$s
      exposure time, and the 2 arcsec slit.

* FIT_TEST_INTENSITIES_FE_13: This routine finds the best-fit density and path length for a given
  set of fake intensities. For example,
	
```
       model log_n = 9.76 +- 0.013  [9.76]
      model log_ds = 7.95 +- 0.027  [7.95]
              chi2 = 2.8
           Line    Imodel      Iobs    SigmaI  dI/Sigma      dI/I
        196.525    463.45    462.16      7.41      0.17       0.3
        200.021    500.55    508.83     10.08      0.82       1.6
        201.121    555.42    550.57     12.81      0.38       0.9
        202.044    796.81    787.23     18.37      0.52       1.2
        203.165    216.98    230.17     11.83      1.12       5.7
        203.826   2511.51   2504.69     44.59      0.15       0.3
        209.916    137.61    149.46     18.76      0.63       7.9
```

 Looping over all of the realizations of CHIANTI yields a distribution that looks like this.

![The distribution of inferred density and path lengths for a single set of input intensities. Note
 that 1-$\sigma$ errors have been added to the intensities.](fit_test_intensities_fe_13.jpg)
