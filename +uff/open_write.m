function open_write(filename)
    attr_details.Name = 'version';
    attr_details.AttachedTo = '/';
    attr_details.AttachType = 'group';
    hdf5write(filename, attr_details, 'v1.0.0');
end

