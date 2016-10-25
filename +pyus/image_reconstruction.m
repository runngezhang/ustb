function elapsed_time = image_reconstruction(file, recons, implem, method)
data_path = which(file);

ustbdata = py.pyus.Interop.USTBInterop.USTBInteropFilter();
ustbdata.param.update(pyargs('path',data_path));

pyusrecons = pyus.create_pyusrecons(recons);

tic;
bf = py.pyus.Filters.Beamforming.BeamformingFilters.BeamformingFilter();
bf.set_inputs(ustbdata);
bf.param.update(pyargs('recons', pyusrecons, 'method', method));

recons.data = pyus.fromNdarray(bf.data);
elapsed_time=toc;

recons.format = E.signal_format.(char(ustbdata.attr{'format'}));
end

