
#### ISSIAtomicData/phase2_20161006/03_MHD

* mhd2ints: Reads in a PLUTO snapshot, computes Fe XIII intensities for each slice in z, sums the
  results along z to create an integrated intensity.  A small region is selected and the peak and mean intensities are determined. Some nuances

    - The emissivities computed 01_chianti_errors need some 'extra' factors (e.g., elemental
      abundance, ionization fraction) to be useful for comparison with observations in absolute
      units. They are added in here.

![The integrated (in z) intensities in two Fe XIII lines of interest.](mhd2ints.1.jpg)

![The temperature and density along z for the peak intensity in Fe XIII 203.826.](mhd2ints.2.jpg)

The peak and mean Fe XIII 203.826 intensities for the small field of view are

```
peak =    109.1   /   141.0   =    0.77
mean =     59.7   /    95.3   =    0.63
```

* To do: Compute intensities for all of the Fe XIII emission lines of interest. Only the 202 and 204 lines have been calculated at present.
