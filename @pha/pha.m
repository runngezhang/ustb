classdef pha < us_dataset
    %PHA Phased array dataset
    %
    % Dataset to do classical phased array imaging
    
    properties
        angle  %Transmit focus beam angles 
    end
    
    methods
    end
    
    %% HUFF
    methods (Access = public)
        function huff_write(h,filename,group_name)
            %HUFF_DUMP    Dumps all the information of the ultrasound
            %dataset to a group in a HDF5
            %
            %   Syntax:
            %   HUFF_DUMP(file_name,group_name)
            %       file_name                    Name of the hdf5 file
            %       group_name                   Name of the group destination
            %
            %   See also STA, US_DATASET
            
            % call superclass function
            h.huff_write@us_dataset(filename,group_name);
            
            % dataset subtype
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
                H5T.enum_insert (filetype, 'STA', 0); 
                H5T.enum_insert (filetype, 'CPW', 1); 
                H5T.enum_insert (filetype, 'VS', 2); 
                H5T.enum_insert (filetype, 'BS', 3); 
                H5T.enum_insert (filetype, 'PHA', 4); 
            gid = H5G.open(file,group_name);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'subtype', filetype, space, 'H5P_DEFAULT');
            
            H5A.write (attr, filetype, uint32(4)); % <---- TYPE PHA
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);
            
            % add angles
            dset_details.Location = group_name;
            dset_details.Name = 'angle';
            hdf5write(filename, dset_details, h.angle, 'WriteMode', 'append');
        end
        function huff_read(h,filename,location)
            %HUFF_READ    Dumps all the information of the ultrasound
            %dataset to a group in a HDF5
            %
            %   Syntax:
            %   HUFF_READ(file_name,location)
            %       file_name                    Name of the hdf5 file
            %       location                   Name of the group destination
            %
            %   See also US_DATASET
            
            % call superclass function
            h.huff_read@us_dataset(filename,location);
            
            % subtype
            dataset_subtype=h5readatt(filename,location,'subtype');
            assert(strcmp(dataset_subtype,'PHA'),'Ultrasound signal subtype does not match!');
            
            % read angles
            h.angle=h5read(filename,[location '/angle']);
        end
    end
    
end

