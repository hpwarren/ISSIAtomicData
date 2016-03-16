
pro nrl_save_hdf_variable, var, name, file_id, verbose=verbose

  if n_elements(var) eq 1 then var = [var]
  dtype_id = h5t_idl_create(var)
  dspace_id = h5s_create_simple(size(var,/dimensions))
  dset_id = h5d_create(file_id, name, dtype_id, dspace_id)
  h5d_write, dset_id, var

  class = H5T_GET_CLASS(H5D_GET_TYPE(dset_id))
  if keyword_set(verbose) then print,' saving variable : '+name+' as '+class

  h5s_close, dspace_id
  h5d_close, dset_id

end

pro nrl_save_hdf_string, var, name, file_id, verbose=verbose

  len = strlen(var)
  u = uniq(len,sort(len))
  if keyword_set(verbose) then print,'  nUnique strlen = '+trim(n_elements(u))
  if n_elements(u) eq 1 then begin
    if keyword_set(verbose) then print,'  fixed length string'
    if n_elements(var) eq 1 then var = [var]
    dtype_id = h5t_idl_create(var)
    dspace_id = h5s_create_simple(size(var,/dimensions))
    dset_id = h5d_create(file_id, name, dtype_id, dspace_id)
    h5d_write, dset_id, var
  endif else begin
    if keyword_set(verbose) then print,'  variable length string, saving to vlen'
    vArray = h5t_str_to_vlen(var)
    dtype_id = h5t_vlen_create(h5t_idl_create(var[0]))
    dspace_id = h5s_create_simple(n_elements(vArray))
    dset_id = h5d_create(file_id, name, dtype_id, dspace_id)
    h5d_write, dset_id, vArray
  endelse

  class = H5T_GET_CLASS(H5D_GET_TYPE(dset_id))
  if keyword_set(verbose) then print,' saving string : '+name+' as '+class

  h5s_close, dspace_id
  h5d_close, dset_id

end

pro nrl_save_hdf_structure, var, name, file_id, verbose=verbose

  dtype_id = h5t_idl_create(var)
  dspace_id = h5s_create_simple(n_elements(var))
  dset_id = h5d_create(file_id, name, dtype_id, dspace_id)
  h5d_write, dset_id, var

  class = H5T_GET_CLASS(H5D_GET_TYPE(dset_id))
  if keyword_set(verbose) then print,' saving structure : '+name+' as '+class

  h5s_close, dspace_id
  h5d_close, dset_id

end

pro nrl_save_hdf, _EXTRA=e, file=file, text=text, verbose=verbose

  if n_elements(e) eq 0 then return

  var_names = strlowcase(tag_names(e))

  ;; --- deal with file name 

  if keyword_set(file) then fileout = file else fileout = 'save'
  ext = '.h5'
  fileout = file_dirname(fileout) + path_sep() + $
      file_basename(fileout, ext, /fold_case) + ext
  if file_test(fileout) then file_delete, fileout

  file_id = h5f_create(fileout)
  if keyword_set(verbose) then print,' saving to file = '+fileout
  
  ;; --- deal with text

  if keyword_set(text) then textin = text else textin = 'no user supplied text'
  if n_elements(textin) gt 1 then textin = strjoin(textin, 10b, /single)
  nrl_save_hdf_string, textin, 'text', file_id

  ;; --- deal with variables

  foreach name, var_names, i do begin
    var = e.(i)
    type = datatype(var)
    case type of
      'STR': nrl_save_hdf_string, var, name, file_id, verbose=verbose
      'STC': nrl_save_hdf_structure, var, name, file_id, verbose=verbose
       else: nrl_save_hdf_variable, var, name, file_id, verbose=verbose
    endcase
  endforeach

  ;; --- save a list of variable names for later use
  var_names = strjoin(var_names, ',', /single)
  nrl_save_hdf_string, var_names, 'nrl_save_var_names', file_id

  h5f_close, file_id

return
end
