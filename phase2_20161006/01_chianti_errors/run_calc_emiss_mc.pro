
pro run_calc_emiss_mc

  tic

  sngl_ion = 'fe_13'
  wvl_list = [196.525, 200.021, 201.121, 202.044, 203.165, 203.826, 209.916]
  nsim = 1000
  normal = 1
  out_file = 'fe_13.monte_carlo'
  dl = 0.1

  calc_emiss_mc, sngl_ion, wvl_list, dl=dl, nsim=nsim, normal=normal, $
                 out_file=out_file

  toc

end
