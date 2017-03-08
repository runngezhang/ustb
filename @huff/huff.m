classdef huff < handle
%HUFF    Class for HUFF reading and writting.
%
%   See also HUFF/HUFF, US_DATASET

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/02/03 $
    
    properties (SetAccess = public)
        filename                % name of the filename
        permission='r'          % String definiting the type of operation sought. Use 'r' for read-only, 'w' for creating or overwritting, and 'rw' for adding  information to an existing file. (Default = 'r')
        ultrasound_signal       % list of ultrasound signal datasets
        spatial_reconstruction  % list of spatial reconstruction dataset
    end
    
    properties (SetAccess = protected)
        version='0.0.2'         % version of the HUFF specification
    end
    
    %% constructor
    methods (Access = public)        
        function h = huff(filename,permission)
            %HUFF    Constructor of the HUFF class.
            %
            %   Syntax:
            %   HUFF(filename) 
            %       filename         Name of the file to be read/written
            %       permission       String definiting the type of
            %                        operation sought. Use 'r' for
            %                        read-only, 'w' for creating or
            %                        overwritting, and 'rw' for adding 
            %                        information to an existing file. 
            %                        (Default = 'r')
            %
            %   See also HUFF
              
            h.filename = filename;
            if exist('permission') h.permission=permission; end
            
            % we check if the operation is feasible and create the file in
            % case it is not there
            switch(h.permission)
                case 'w'
                    fcpl = H5P.create('H5P_FILE_CREATE');
                    fapl = H5P.create('H5P_FILE_ACCESS');
                    fid = H5F.create(filename,'H5F_ACC_TRUNC',fcpl,fapl);
                    H5F.close(fid);
                case 'r'
                    fid = H5F.open(filename);
                    H5F.close(fid);
                case 'rw'
                    fid = H5F.open(filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                    H5F.close(fid);
                otherwise
                    error('Unknown permission!');
            end
        end
    end
    
    %% append
    methods (Access = public)
        function stage(h,dataset)
            %STAGE    Method that copies the handle of a known dataset into a list to add to a hdf5 file
            %
            %   Syntax:
            %   STAGE(dataset)
            %       dataset          Instance of a known dataset: (sta, cpw, vs, or reconstruction)
            %
            %   See also HUFF
            known_datasets = {'sta', 'cpw', 'vs','pha'};
            if ismember(class(dataset),known_datasets)
                n=length(h.ultrasound_signal);
                h.ultrasound_signal{n+1}=dataset;
            elseif strcmp(class(dataset),'reconstruction') 
                n=length(h.spatial_reconstruction);
                h.spatial_reconstruction{n+1}=dataset;
            else
                error('Unknown dataset!');
            end
        end
    end
    
    %% write
    methods (Access = public)
        function write(h)
            %WRITE    Method that writes the listed datasets into the specified hd5f
            %
            %   Syntax:
            %   WRITE()
            %
            %   See also HUFF
            
            
            % write version
            attr_details.Name = 'version';
            attr_details.AttachedTo = '/';
            attr_details.AttachType = 'group';
            hdf5write(h.filename, attr_details, h.version);
            
            % We create the /US metagroup in case it is not there
            try
                h5info(h.filename,'/US')
            catch 
                fid = H5F.open(h.filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                gid = H5G.create(fid,'US','H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
                H5G.close(gid);
                H5F.close(fid);
            end
            
            % ultrasound signal dataset loop
            wb=waitbar(0,'Writting ultrasound signal datasets');
            for n=1:length(h.ultrasound_signal)
                waitbar(n/length(h.ultrasound_signal));
                % create new group in metagroup
                group_name=sprintf('US%d',n);
                fid = H5F.open(h.filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                gid = H5G.open(fid,'/US');
                s_gid = H5G.create(gid,group_name,'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
                H5G.close(s_gid);
                H5G.close(gid);
                H5F.close(fid);

                % dump 
                h.ultrasound_signal{n}.huff_write(h.filename,['US/' group_name]);
            end
            close(wb);
            
            % We create the /SR metagroup in case it is not there
            try
                h5info(h.filename,'/SR')
            catch 
                fid = H5F.open(h.filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                gid = H5G.create(fid,'SR','H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
                H5G.close(gid);
                H5F.close(fid);
            end
            
            % image reconstruction dataset loop
            wb=waitbar(0,'Writting spatial reconstruction datasets');
            for n=1:length(h.spatial_reconstruction)
                waitbar(n/length(h.spatial_reconstruction));
                % create new group in metagroup
                group_name=sprintf('SR%d',n);
                fid = H5F.open(h.filename,'H5F_ACC_RDWR','H5P_DEFAULT');
                gid = H5G.open(fid,'/SR');
                s_gid = H5G.create(gid,group_name,'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
                H5G.close(s_gid);
                H5G.close(gid);
                H5F.close(fid);

                % dump 
                h.spatial_reconstruction{n}.huff_write(h.filename,['SR/' group_name]);
            end
            close(wb);
        end
    end
    
    %% write
    methods (Access = public)
        function read(h)
            %READ    Method that reads all the datasets in the hdf5 and creates an object of the corresponding class in USTB 
            %
            %   Syntax:
            %   READ()
            %
            %   See also HUFF
            
            % Reading ultrasound signal datasets from file
            info=h5info(h.filename,'/US');
            wb=waitbar(0,'Reading ultrasound signal datasets');
            for n=1:length(info.Groups)
                waitbar(n/length(info.Groups));
                location=info.Groups(n).Name;
                dstype=h5readatt(h.filename,location,'type');
                if strcmp(dstype,'US')
                    subtype=h5readatt(h.filename,location,'subtype');
                    switch(subtype{1})
                        case 'STA'
                            h.ultrasound_signal{n}=sta();
                            h.ultrasound_signal{n}.huff_read(h.filename,location);
                        case 'CPW'
                            h.ultrasound_signal{n}=cpw();
                            h.ultrasound_signal{n}.huff_read(h.filename,location);
                        case 'VS'
                            h.ultrasound_signal{n}=vs();
                            h.ultrasound_signal{n}.huff_read(h.filename,location);
                        case 'PHA'
                            h.ultrasound_signal{n}=pha();
                            h.ultrasound_signal{n}.huff_read(h.filename,location);
                        case 'BS'
                            disp('BS is not supported yet');                
                        otherwise
                            warning(sprintf('Unknown subtype %s of ultrasound signal dataset in %s\n',subtype,location));
                    end
                else
                    warning(sprintf('Unknown type %s of ultrasound signal dataset in %s\n',dstype,location));
                end
            end
            close(wb);

            % Reading spatial reconstruction datasets from file
            info=h5info(h.filename,'/SR');
            wb=waitbar(0,'Reading spatial reconstruction datasets');
            for n=1:length(info.Groups)
                waitbar(n/length(info.Groups));
                location=info.Groups(n).Name;
                dstype=h5readatt(h.filename,location,'type');
                if strcmp(dstype,'SR')
                    h.spatial_reconstruction{n}=reconstruction();
                    h.spatial_reconstruction{n}.huff_read(h.filename,location);
                else
                    warning(sprintf('Unknown type %s of ultrasound signal dataset in %s\n',dstype,location));
                end
            end
            close(wb);
        end
    end
end

