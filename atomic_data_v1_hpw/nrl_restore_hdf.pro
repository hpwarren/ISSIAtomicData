
pro nrl_restore_hdf, _REF_EXTRA=e, file=file, text=text, inquire=inquire

  if n_elements(file) eq 0 then file = 'save.h5'
  if file_test(file) eq 0 then begin
    message,' input file not found '+file,/info
    return
  end
  file_id = h5f_open(file)

  dset_id = h5d_open(file_id,'nrl_save_var_names')
  var_names = h5d_read(dset_id)
  h5d_close, dset_id
  var_names = strsplit(var_names,',',/extract)

  if keyword_set(inquire) then begin
    print,' available variables '
    foreach v, var_names do print,'  '+v
    return
  endif

  dset_id = h5d_open(file_id,'text')
  text = h5d_read(dset_id)
  h5d_close, dset_id
  text = strsplit(text, 10b, /extract)

  if n_elements(e) eq 0 then return
  requested_var_names = strlowcase(e)

  foreach name,requested_var_names, i do begin
    match = where(name eq var_names,nMatch)
    if nMatch eq 0 then begin
      message,' requested variable not found : '+name,/info
      nrl_restore_hdf, file=file, /inquire
      return
    endif else begin
      dset_id = h5d_open(file_id,name)
      data = h5d_read(dset_id)
      type = H5T_GET_CLASS(H5D_GET_TYPE(dset_id))
      h5d_close, dset_id
      ;; --- look for vlen arrays
      if type eq 'H5T_VLEN' then data = h5t_vlen_to_str(data)
      if n_elements(data) eq 1 then data = data[0]
      (scope_varfetch(name, /ref_extra)) = data
    endelse
  endforeach

  h5f_close, file_id

return
end
