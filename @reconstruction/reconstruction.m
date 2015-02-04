classdef reconstruction < handle
%RECONSTRUCTION   Reconstruction definition
%
%   See also BEAM, SCAN
    properties  (SetAccess = public)
        name=''                     % String containing the name of the reconstruction
        creation_date=''            % String containing the date the reconstruction was created
        transmit_beam=beam()        % BEAM class defining transmit beam
        receive_beam=beam()         % BEAM class defining receive beam            
        scan                        % SCAN class defining the scan area 
        format=E.signal_format.RF   % format of the signal
        data                        % matrix containing the reconstructed raw signal
        envelope                    % matrix containing the envelope of the reconstructed signal
        frames                      % Number of frames in the reconstruction
        central_frequency           % central frequency [Hz]
        bandwidth                   % central frequency [Hz]
    end
    
    %% constructor
    methods (Access = public)
        function h = reconstruction(name, object)
            %RECONSTRUCTION    Constructor of the RECONSTRUCTION class.
            %
            %   Syntax:
            %   RECONSTRUCTION(name, object) 
            %       name                    Name of the reconstruction
            %       object                  Instance of a reconstruction class to copy values from
            %
            %   See also BEAM, LINEAR_SCAN

            t_now=now;
            h.creation_date=[date sprintf('-%d-%d-%d',hour(t_now),minute(t_now),round(second(t_now)))];

            if exist('object') 
                h.copy(object);
            end

            if exist('name') 
                h.name=name;
            end
        end
    end
    
    %% copy
    methods (Access = public)
        function copy(h,object)
            %COPY    Copy the values from another reconstruction
            %
            %   Syntax:
            %   COPY(object) 
            %       object       Instance of a reconstruction class
            %
            %   See also BEAM, LINEAR_SCAN
            assert(isa(object,class(h)),'Class of the input object is not identical'); 

            h.name=object.name;     
            h.transmit_beam=object.transmit_beam;   % beam is not a handle class thus we can copy it
            h.receive_beam=object.receive_beam;     % beam is not a handle class thus we can copy it       
            h.scan=object.scan;                     % scan is not a handle class thus we can copy it       
            
            h.format=object.format;                 
            h.data=object.data;
            h.envelope=object.envelope;
            h.frames=object.frames;
            h.central_frequency=object.central_frequency;     
            h.bandwidth=object.bandwidth;
        end
    end
   
    %% calculate apodization
    methods (Access = public)
        function [apo]= calculate_apodization(h,beam,x0)
            %CALCULATE_APODIZATION    Calculates apodization for a given set of elements
            %
            %   Syntax:
            %   calculate_apodization(beam, elements_x_coordinate)
            %       BEAM                    Class specifying the beam
            %       elements_x_coordinate   Vector with the x coordinates of the elements
            %
            %   See also BEAM, LINEAR_SCAN
            
            % checking format
            assert(size(x0,1)==h.scan.pixels, 'The element position vector must be a column vector!');
            
            % we take the number of apodization values from x0
            number_transmitting_events=size(x0,2);
            
            % apodization computation 
            if beam.apodization==E.apodization_type.none
                apo=ones(h.scan.pixels,number_transmitting_events);
            else
                apo=zeros(h.scan.pixels,number_transmitting_events);
                Aperture=h.scan.z./beam.f_number;
                for n=1:number_transmitting_events
                    %xd=abs(x0(:,n)-h.scan.x+h.scan.z*tan(beam.steer_angle));
                    xd=h.scan.lateral_distance(x0(:,n),beam.steer_angle);
                    switch(beam.apodization)
                        case E.apodization_type.boxcar
                            apo(:,n)=double(xd<Aperture/2); 
                        case E.apodization_type.hanning
                            apo(:,n)=double(xd<Aperture/2).*(0.5 + 0.5*cos(2*pi*xd./Aperture)); 
                        case E.apodization_type.tukey25
                            roll=0.25;
                            apo(:,n)=(xd<(Aperture/2*(1-roll))) + (xd>(Aperture/2*(1-roll))).*(xd<(Aperture/2)).*0.5.*(1+cos(2*pi/roll*(xd./Aperture-roll/2-1/2)));                               
                        case E.apodization_type.tukey50
                            roll=0.5;
                            apo(:,n)=(xd<(Aperture/2*(1-roll))) + (xd>(Aperture/2*(1-roll))).*(xd<(Aperture/2)).*0.5.*(1+cos(2*pi/roll*(xd./Aperture-roll/2-1/2)));                               
                        case E.apodization_type.tukey80
                            roll=0.8;
                            apo(:,n)=(xd<(Aperture/2*(1-roll))) + (xd>(Aperture/2*(1-roll))).*(xd<(Aperture/2)).*0.5.*(1+cos(2*pi/roll*(xd./Aperture-roll/2-1/2)));                               
                        otherwise
                            error('Unknown apodization type!');
                    end
                end
            end
            
            % our implementation of edge smoothing -> inspired by tukey window
            if(size(apo,2)>(2*(beam.smoothing+1)))
                % compute edge smoothing mask
                mask=0.5-0.5*cos((0:beam.smoothing)/beam.smoothing*pi);  % vector mask
                mask=(ones(h.scan.pixels,1)*mask(2:(beam.smoothing+1))); % matrix mask
                
                % modifying apodization edges
                apo(:,1:beam.smoothing)=apo(:,1:beam.smoothing).*mask;
                apo(:,(number_transmitting_events-beam.smoothing+1):number_transmitting_events)=apo(:,(number_transmitting_events-beam.smoothing+1):number_transmitting_events).*fliplr(mask);
            end
        end
    end
    
    %% HUFF
    methods (Access = public)
        function huff_write(h,filename,location)
            %HUFF_WRITE    Dumps all the information of the ultrasound
            %dataset to a group in a HDF5
            %
            %   Syntax:
            %   HUFF_WRITE(file_name,location)
            %       file_name                    Name of the hdf5 file
            %       location                   Name of the group destination
            %
            %   See also RECONSTRUCTION
            
            % Dataset type
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
                H5T.enum_insert (filetype, 'US', 0); 
                H5T.enum_insert (filetype, 'SR', 1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'type', filetype, space, 'H5P_DEFAULT');
            H5A.write (attr, filetype, uint32(1));  % <--- SR
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);
                        
            % Signal format 
            file = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            filetype = H5T.enum_create ('H5T_NATIVE_INT');
                H5T.enum_insert (filetype, 'RF', 0); 
                H5T.enum_insert (filetype, 'IQ', 1); 
            gid = H5G.open(file,location);
            space = H5S.create_simple (1,1,[]);

            attr = H5A.create (gid, 'signal_format', filetype, space, 'H5P_DEFAULT');
            switch(h.format)
                case E.signal_format.RF
                    H5A.write (attr, filetype, uint32(0));
                case E.signal_format.IQ
                    H5A.write (attr, filetype, uint32(1));
            end
            H5A.close (attr);
            H5G.close(gid);    
            H5S.close (space);
            H5T.close (filetype);
            H5F.close (file);

            % add name
            attr = h.name;
            attr_details.Name = 'name';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');

            % add date
            attr = h.creation_date;
            attr_details.Name = 'creation_date';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');

            % add frames
            attr = h.frames;
            attr_details.Name = 'frames';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');
            
            % add central frequency
            attr = h.central_frequency;
            attr_details.Name = 'central_frequency';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');
            
            % add bandwidth
            attr = h.bandwidth;
            attr_details.Name = 'bandwidth';
            attr_details.AttachedTo = location;
            attr_details.AttachType = 'group';
            hdf5write(filename, attr_details, attr, 'WriteMode', 'append');
            
            % add data
            dset_details.Location = [location '/data'];
            dset_details.Name = 'real';
            hdf5write(filename, dset_details, real(h.data), 'WriteMode', 'append');
            dset_details.Name = 'imag';
            hdf5write(filename, dset_details, imag(h.data), 'WriteMode', 'append');

            % add envelope
            dset_details.Location = [location];
            dset_details.Name = 'envelope';
            hdf5write(filename, dset_details, h.envelope, 'WriteMode', 'append');
            
            % create new groups for substructures
            fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
            gid = H5G.open(fid,location);
            gid2 = H5G.create(gid,'transmit_beam','H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
            gid3 = H5G.create(gid,'receive_beam','H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
            gid4 = H5G.create(gid,'scan','H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
            H5G.close(gid4);
            H5G.close(gid3);
            H5G.close(gid2);
            H5G.close(gid);
            H5F.close(fid);

            % dump substructures
            h.transmit_beam.huff_write(filename,[location '/transmit_beam']);
            h.receive_beam.huff_write(filename,[location '/receive_beam']);
            h.scan.huff_write(filename,[location '/scan']);

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
            
            % type
            dataset_type=h5readatt(filename,location,'type');
            assert(strcmp(dataset_type,'SR'),'Spatial reconstruction type does not match!');
            
            % read signal format 
            signal_format=h5readatt(filename,location,'signal_format');
            switch(signal_format{1})
                case 'RF'
                    h.format=E.signal_format.RF;
                case 'IQ'
                    h.format=E.signal_format.IQ;
                otherwise
                    error('Unknown signal format!');
            end

            % read substructures
            h.transmit_beam=h.transmit_beam.huff_read(filename,[location '/transmit_beam']);
            h.receive_beam=h.receive_beam.huff_read(filename,[location '/transmit_beam']);
            scan_type=h5readatt(filename,[location '/scan'],'scan_type');
            switch (scan_type{1})
                case 'linear_scan'
                    h.scan=linear_scan();
                    h.scan=h.scan.huff_read(filename,[location '/scan']);  
                case 'sector_scan'
                    h.scan=sector_scan();
                    h.scan=h.scan.huff_read(filename,[location '/scan']);  
                otherwise
                    warning('Unknown scan format!');
            end
            
            % read frames
            h.frames=h5readatt(filename,location,'frames');

            % read creation_date
            h.creation_date=h5readatt(filename,location,'creation_date');
            
            % read central_frequency
            h.central_frequency=h5readatt(filename,location,'central_frequency');

            % read bandwidth
            h.bandwidth=h5readatt(filename,location,'bandwidth');

            % read data
            real_part=h5read(filename,[location '/data/real']);
            imag_part=h5read(filename,[location '/data/imag']);
            h.data=real_part+1i*imag_part;
            
            % read envelope --> optional
            h.envelope=h5read(filename,[location '/envelope']);
            
        end
    end
    
    %% presentation methods
    methods (Access = public)    
        function im=calculate_envelope(h)
            %CALCULATE_ENVELOPE    Calculates the envelope of the beamformed signal
            %
            %   Syntax:
            %   envelope_signal=calculate_envelope(data)
            %       envelope_signal         Matrix with the envelope with the same dimensions of data
            %
            %   See also RECONSTRUCTION
            
            switch(h.format)
                case E.signal_format.RF
                    dz=h.scan.depth_step();
                    Fs=1540/2/dz;                   % effective sampling frequency (Hz)
                    maximum_frequency=h.central_frequency+h.bandwidth;
                    ratio=Fs/maximum_frequency;
                    
                    im=zeros(size(h.data));
                    
                    if(ratio<4)
                        warning(sprintf('The spatial sampling in the beam direction does not allow to compute the Hilbert transform on pixel data. Displayed envelopes will be innacurate. To solve this issue please consider: working with demodulated signal (IQ), or increasing the sampling resolution in the beam direction by a factor of %0.2f',4/ratio));
                    end
                    for f=1:h.frames
                        im(:,:,f)=abs(hilbert(h.data(:,:,f)));
                    end
                    
                case E.signal_format.IQ
                    im=abs(h.data);
                otherwise
                    error('Unknown signal format!');
            end
        end
        
        function im=show(h,dynamic_range)
            %SHOW    Plots the envelope of the beamformed data and returns a copy of the image
            %
            %   Syntax:
            %   image=show(dynamic_range)
            %       image           Matrix with the envelope of the normalised bemformed signal on dB
            %       dynamic_range   Desired dynamic range of the displayed images
            %
            %   See also RECONSTRUCTION
       
            if ~exist('dynamic_range') dynamic_range=60; end
            
            % computing envelope
            if isempty(h.envelope) h.envelope=h.calculate_envelope(); end
            
            % Ploting image reconstruction
            im=20*log10(h.envelope./max(h.envelope(:)));
            figure; set(gca,'fontsize',16); 
            for f=1:h.frames
                x_lim=[min(h.scan.x_matrix(:)) max(h.scan.x_matrix(:))]*1e3;
                z_lim=[min(h.scan.z_matrix(:)) max(h.scan.z_matrix(:))]*1e3;
                % black background
                %pcolor(x_lim,z_lim,[-dynamic_range -dynamic_range; -dynamic_range -dynamic_range]); shading flat; colormap gray; caxis([-dynamic_range 0]); colorbar; hold on;
                pcolor(h.scan.x_matrix*1e3,h.scan.z_matrix*1e3,im(:,:,f)); shading flat; colormap gray; caxis([-dynamic_range 0]); colorbar; hold on;
                axis equal manual;
                xlabel('x [mm]');
                ylabel('z [mm]');
                set(gca,'YDir','reverse');
                set(gca,'fontsize',16);
                axis([x_lim z_lim]);
                title(sprintf('%s (%s)',char(h.name),char(h.format))); 
                drawnow; hold off;
                pause(0.05);
            end
        end
        
        function write_video(h,filename,dynamic_range,figure_size)
            %MAKE_VIDEO    Creates and stores a video with the beamformeddataset
            %
            %   Syntax:
            %   write_video(filename,dynamic_range)
            %       filename        Path and name where the video will be written
            %       dynamic_range   Desired dynamic range of the displayed images (optional)
            %       figure_size     Size of the exported video in pixels [horizontal, vertical] (optional)
            %
            %   See also RECONSTRUCTION
       
            if ~exist('dynamic_range') dynamic_range=60; end
            if ~exist('figure_size') figure_size=[1024 768]; end
            
            
            %% making a video -> envelope
            writerObj = VideoWriter(filename);
            open(writerObj);
            
            % computing envelope
            if isempty(h.envelope) h.envelope=h.calculate_envelope(); end
            
            im=20*log10(h.envelope./max(h.envelope(:)));
            figure; set(gca,'fontsize',16); 
            for f=1:h.frames
                x_lim=[min(h.scan.x_matrix(:)) max(h.scan.x_matrix(:))]*1e3;
                z_lim=[min(h.scan.z_matrix(:)) max(h.scan.z_matrix(:))]*1e3;
                pcolor(h.scan.x_matrix*1e3,h.scan.z_matrix*1e3,im(:,:,f)); shading flat; colormap gray; caxis([-dynamic_range 0]);colorbar; hold on;
                axis equal manual;
                xlabel('x [mm]');
                ylabel('z [mm]');
                set(gca,'YDir','reverse');
                set(gca,'fontsize',16);
                axis([x_lim z_lim]);
                title(sprintf('Frame %d',f)); 
                set(gcf,'Position',[200   100   figure_size(1) figure_size(2)]);
                drawnow;

                frame = getframe(gcf,[0 0 figure_size(1) figure_size(2)]);
                writeVideo(writerObj,frame);

                hold off;
            end

            close(writerObj);
            
        end
    end

    %% Set methods
    methods
        function h=set.data(h,input_data)
            % column-based format
            if(size(input_data,1)==h.scan.pixels)
                h.frames=size(input_data,2);
            
                % reshape beamformed image
                h.data=reshape(input_data,[size(h.scan.x_matrix,1) size(h.scan.x_matrix,2) h.frames]);
                h.envelope=[];
            % matrix format
            else
                % check that the data format matches the specification
                assert((size(input_data,1)==size(h.scan.x_matrix,1))&&(size(input_data,2)==size(h.scan.x_matrix,2)), 'The number of rows in the data does not match the scan specification!');
                h.frames=size(input_data,3);
            
                % copy beamformed image
                h.data=input_data;
                h.envelope=[];
            end
        end
    end

end

