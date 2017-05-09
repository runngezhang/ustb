function append(filename, location, name, value)
    if isa(value,'double')
        if isreal(value)
            dset_details.Location = location;
            dset_details.Name = name;
            hdf5write(filename, dset_details, single(value), 'WriteMode', 'append');
        else
            dset_details.Location = [location '/' name];
            dset_details.Name = 'real';
            hdf5write(filename, dset_details, single(real(value)), 'WriteMode', 'append');
            dset_details.Name = 'imag';
            hdf5write(filename, dset_details, single(imag(value)), 'WriteMode', 'append');
        end
    else
        error('Class not supported yet');
    end       
end

