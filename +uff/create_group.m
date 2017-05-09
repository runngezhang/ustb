function new_location=create_group(filename, location, group_name, version)
    group_name = 'channel_data';
    fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
    if isempty(location)
        gid = H5G.open(fid,'/');
    else
        gid = H5G.open(fid,location);
    end
    s_gid = H5G.create(gid,group_name,'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
    H5G.close(s_gid);
    H5G.close(gid);
    H5F.close(fid);
    new_location = [location '/' group_name];
            
   % write group version 
   attr_details.Name = 'version';
   attr_details.AttachedTo = new_location;
   attr_details.AttachType = 'group';
   hdf5write(filename,attr_details,version,'WriteMode','append');
end

