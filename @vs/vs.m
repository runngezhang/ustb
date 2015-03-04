classdef vs < us_dataset
%VS    Virtual source (vs) dataset.
%
%   See also US_DATASET, VS.VS

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/01/29 $
    
    properties (SetAccess = public)
        source       % vector containing the angles (rad)
    end
    
    %% constructor
    methods (Access = public)        
        function h = vs(name,input_format,input_c0,input_source,input_time,input_data,input_geom,input_modulation_frequency)
            %VS    Constructor of the VS class.
            %
            %   Syntax:
            %   VS(name,format,c0,source,time,data,geom,modulation_frequency) 
            %       name                    Name of the dataset
            %       format                  Format of the signal (E.signal_format.RF, default=E.signal_format.IQ)
            %       c0                      Reference speed of sound (m/s)
            %       source                  Location of virtual sources [x, y, z] (m)
            %       time                    Time vector (s)
            %       data                    Numerical data [time_samples, channels, firings, frames]
            %       geom                    Probe geometry [x, y, z] (m)
            %       modulation_frequency    Modulation frequency (Hz) - required only for IQ format
            %
            %   See also VS, US_DATASET
              
            % we allow object instantiation without parameters
            h@us_dataset(); 
            
            if exist('name') 
                % required data
                h.name = name;
                h.format = input_format;
                h.c0 = input_c0;
                h.time = input_time;
                h.data = input_data;
                h.geom = input_geom;
                h.source = input_source;

                if(h.format==E.signal_format.IQ) 
                    h.modulation_frequency=input_modulation_frequency;
                end

                % checks
                assert(h.channels==size(h.geom,1),'The length of the geometry vector should match number of channels in the dataset');
                assert(h.firings==size(h.source,1),'The source definition should match the number of firings in the dataset');
            end
        end

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
            gid = H5G.open(file,group_name);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'subtype', filetype, space, 'H5P_DEFAULT');
            
            H5A.write (attr, filetype, uint32(2)); % <---- TYPE VS

            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);
            
            % add transducer geometry
            dset_details.Location = group_name;
            dset_details.Name = 'source';
            hdf5write(filename, dset_details, h.source, 'WriteMode', 'append');
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
            assert(strcmp(dataset_subtype,'VS'),'Ultrasound signal subtype does not match!');
            
            % read angles
            h.source=h5read(filename,[location '/source']);
        end
    end

    %% image reconstruction
    methods (Access = public)        
        function elapsed_time = image_reconstruction(h,recons,implem)
            %IMAGE_RECONSTRUCTION    Method for image reconstruction of virtual source datasets
            %
            %   Syntax:
            %   IMAGE_RECONSTRUCTION(transmit_beam,receive_beam,scan,implementation)
            %       reconstruction          Class containing the reconstruction specification
            %       implementation          Enumeration specifying the algorithm 
            %
            %   See also RECONSTRUCTION, VS
             
            % default implementation
            if ~exist('implem') implem=E.implementation.mex; end
            
            % initial time -> performance test
            tic; 
            
            % loop over orientations
            total_data=zeros(size(recons.scan.x_matrix,1),size(recons.scan.x_matrix,2),length(recons.orientation),size(h.data,4));
            temporal_data=zeros(recons.scan.pixels,size(h.data,4));
            for o=1:length(recons.orientation)
                % precompute transmit and receive apodization
                xT=recons.scan.x*ones(1,h.firings)-(recons.scan.z*ones(1,h.firings)).*(ones(recons.scan.pixels,1)*h.source(:,1).'-recons.scan.x*ones(1,h.firings))./(ones(recons.scan.pixels,1)*h.source(:,3).'-recons.scan.z*ones(1,h.firings)); % position of equivalent receive element -> Alfonso's equation 
                valid_apodization=(xT>h.geom(1,1))&(xT<h.geom(end,1));            % check we don't get out of the aperture
                h.transmit_apodization = valid_apodization.*recons.calculate_apodization(recons.orientation(o).transmit_beam,xT);
                xR=ones(recons.scan.pixels,1)*(h.geom(:,1).');                    % position of receive element
                h.receive_apodization = recons.calculate_apodization(recons.orientation(o).receive_beam,xR);

                % launch selected implementation
                temporal_data=h.launch_implementation(recons,implem);

                % reshape matrix to include orientation dimensions
                total_data(:,:,o,:)=reshape(temporal_data,[size(recons.scan.x_matrix,1) size(recons.scan.x_matrix,2) 1 size(temporal_data,2)]);
            end
            
            % copy data to reconstruction
            recons.data=total_data;
            recons.format=h.format;
            if(h.format==E.signal_format.RF)
                [recons.central_frequency recons.bandwidth]=tools.estimate_frequency(h.time,h.data);
            end
            
            % elapsed time for performance test
            elapsed_time=toc; 
        end
    end
    
    %% set methods
    methods  
        function set.source(h,input_geom)
            assert(ndims(input_geom)==2&&size(input_geom,2)==3, 'Wrong source definition. It should be a three column matrix [x, y, z]')
            h.source=input_geom;
        end
    end
end

