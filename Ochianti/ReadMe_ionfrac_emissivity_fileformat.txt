Emissivities are in the files emissivity_OJ_XXYYYY.txt, where J is the ionization state and XX.YYYY is the wavelength of the transition in Angstroms
	first row: atomic number, ionization state, wavelength, size of temperature array, size of density array
	second row: first element is size of temperature array repeated, followed by the temperature array in log10(T[K])
	subsequent rows: first element is the log10(density[cm^-3]), followed by the bare emissivities without ionization balance multiplied in [ph cm^3/s]

Ionization balances are in the files ionfrac_OJ.txt, where J is the ionization state
	first row: atomic number, ionization state, size of temperature array
	second row: temperature array in log10(T[K])
	third row: ionization fraction

The emissivities are generated using CHIANTI v7.1.3
Ionization fractions are from CHIANTI/dbase/ioneq/chianti.ioneq

Multiply emissivity at a given temperature by the ionization fraction at the same temperature before using it to estimate line flux.
