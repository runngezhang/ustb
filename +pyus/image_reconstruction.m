function [elapsed_time, res_data] = image_reconstruction(file, recons, implem, method)
data_path = which(file);

ustbdata = py.pyus.Interop.USTBInterop.USTBInteropFilter();
ustbdata.param.update(pyargs('path',data_path));
ustbdata.update();

if isequal(class(recons), class(py.pyus.Interop.USTBInterop.Reconstruction))
    pyusrecons = recons;
else
    pyusrecons = pyus.create_pyusrecons(recons);
end

tic;
bf = py.pyus.Filters.Beamforming.BeamformingFilters.BeamformingFilter();
bf.set_inputs(ustbdata);
bf.param.update(pyargs('recons', pyusrecons, 'method', method));

res_data = pyus.from_ndarray(bf.data);

if ~isequal(class(recons), class(py.pyus.Interop.USTBInterop.Reconstruction))
    recons.data = res_data';
    recons.format = E.signal_format.(char(ustbdata.attr{'format'}));
end
elapsed_time=toc;
    
end

