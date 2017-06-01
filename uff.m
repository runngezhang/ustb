classdef uff
    %UFF   defines the read/write functions for Ultrasound File Format (UFF)
    % structures
    %
    %   See also CHANNEL_DATA, BEAMFORMED_DATA
    
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %
    %   $Last updated: 2017/05/10$
    
    % UFF variables
    properties (SetAccess = private)
        filename
        version = 'v1.0.1';
        mode = 'append';
        verbose = true;
    end
    
    %% constructor
    methods (Access = public)
        function h=uff(filename,mode)
            
            if nargin<1
                error('UFF: Missing filename');
            end
            if nargin<2
                mode='append';
            end
            
            switch mode
                case 'write'
                    attr_details.Name = 'version';
                    attr_details.AttachedTo = '/';
                    attr_details.AttachType = 'group';
                    hdf5write(filename, attr_details, h.version);
                case 'read'
                    assert(exist(filename,'file')==2,sprintf('File %s not found. Check the path and filename, or chose ''write'' mode to create it.',filename));

				    % check that versions are the same
                    file_version=h5readatt(filename, '/','version');  % read file version
                    file_version=file_version{1};                       % from cell to string
                    file_version=file_version(int32(file_version)>0);   % removing 0's from 0-terminated strings
                    switch file_version
                        case h.version
                        case 'v1.0.0'
                        otherwise
                            error(sprintf('Unsupported file version (%s). Current UFF version (%s). Please choose a new file instead.',file_version,h.version));
                    end
                case 'append'
                    if(exist(filename,'file')~=2)
                        % create it
                        if h.verbose fprintf('UFF: file %s not found; it shall be created.\n',filename); end
                        attr_details.Name = 'version';
                        attr_details.AttachedTo = '/';
                        attr_details.AttachType = 'group';
                        hdf5write(filename, attr_details, h.version);
                    else
                        % check that versions are the same
                        file_version=h5readatt(filename, '/','version');  % read file version
                        file_version=file_version{1};                       % from cell to string
                        file_version=file_version(int32(file_version)>0);   % removing 0's from 0-terminated strings
                        switch file_version
                            case h.version
                            case 'v1.0.0'
                            otherwise
                                error(sprintf('Unsupported file version (%s). Current UFF version (%s). Please choose a new file instead.',file_version,h.version));
                        end
                    end
                otherwise
                    error(sprintf('Mode %s not supported',mode));
            end
            
            h.mode = mode;
            h.filename = filename;
        end
    end
    
    methods
        function namestr = getobjname(h,obj)
            namestr = inputname(2);
        end
    end
    
    %% Ultrasound File Format
    methods (Access = public)
        function dumped_objects=write(h,object,name,location)
            %% WRITE- Writes object into location
            %
            %   This UFF method writes the object into the specified location
            %   with the provided name. A UFF object must be defined.
            %
            %   dumped_objects=UFF.WRITE(object,name,location)
            %
            %   Example:
            %
            %       channel_data = uff.channel_data();
            %       uff_file = uff('test.uff');
            %       dumped_objects = uff_file.write(channel_data,'chn01');
            %
            %   See also UFF.READ, UFF.UFF, UFF.INDEX
            
            if(nargin<2)||isempty(object)
                error('Missing object to write on UFF file');
            end
            if(nargin<3)||isempty(name)
                error('Missing object name');
            end
            if(nargin<4)
                location=[];
            end
            
            assert(~strcmp(h.mode,'read'),'The file is open in read-only mode');
            
            switch class(object)
                case {'double' 'single'}
                    if isreal(object)
                        h5create(h.filename,[location '/' name], size(object), 'Datatype', 'single', 'ChunkSize',size(object));
                        h5write(h.filename,[location '/' name], single(object));
                        h5writeatt(h.filename,[location '/' name],'class',class(object));
                        h5writeatt(h.filename,[location '/' name],'name',name);
                        h5writeatt(h.filename,[location '/' name],'imaginary',0);
                        h5writeatt(h.filename,[location '/' name],'complex',0);
                        dumped_objects=1;
                    else
                        % real
                        h5create(h.filename,[location '/' name '/real'], size(object), 'Datatype', 'single', 'ChunkSize',size(object));
                        h5write(h.filename,[location '/' name '/real'], single(real(object)));
                        h5writeatt(h.filename,[location '/' name '/real'],'class',class(object));
                        h5writeatt(h.filename,[location '/' name '/real'],'name',name);
                        h5writeatt(h.filename,[location '/' name '/real'],'imaginary',0);
                        
                        % imag
                        h5create(h.filename,[location '/' name '/imag'], size(object), 'Datatype', 'single', 'ChunkSize',size(object));
                        h5write(h.filename,[location '/' name '/imag'], single(imag(object)));
                        h5writeatt(h.filename,[location '/' name '/imag'],'class',class(object));
                        h5writeatt(h.filename,[location '/' name '/imag'],'name',name);
                        h5writeatt(h.filename,[location '/' name '/imag'],'imaginary',1);
                        
                        % group attributes
                        h5writeatt(h.filename,[location '/' name],'class',class(object));
                        h5writeatt(h.filename,[location '/' name],'name',name);
                        h5writeatt(h.filename,[location '/' name],'complex',1);
                        dumped_objects=1;
                    end
                case 'char' 
                    h5create(h.filename,[location '/' name], size(object), 'Datatype', 'single', 'ChunkSize',size(object));
                    h5write(h.filename,[location '/' name], uint16(object));
                    h5writeatt(h.filename,[location '/' name],'class',class(object));
                    h5writeatt(h.filename,[location '/' name],'name',name);
                    dumped_objects=1;
                case 'uff.window'
                    h5create(h.filename,[location '/' name], size(object), 'Datatype', 'single', 'ChunkSize',size(object));
                    h5write(h.filename,[location '/' name], single(object));
                    h5writeatt(h.filename,[location '/' name],'class',class(object));
                    h5writeatt(h.filename,[location '/' name],'name',name);
                    dumped_objects=1;
                case 'cell'
                    % call write for all members in the cell
                    dumped_objects=0;
                    for n=1:numel(object)
                        dumped_objects=dumped_objects+h.write(object{n}, [name '_' sprintf('%04d',n)],[location '/' name]);
                    end
                    
                    % group attributes
                    if dumped_objects
                        h5writeatt(h.filename,[location '/' name],'class',class(object));
                        h5writeatt(h.filename,[location '/' name],'name',name);
                        h5writeatt(h.filename,[location '/' name],'array',1);
                        h5writeatt(h.filename,[location '/' name],'size',size(object));
                    end
                otherwise
                    % UFF structures
                    if (findstr('uff.',class(object)))
                        if numel(object)>1
                            
                            if h.verbose fprintf('UFF: writting %s [%s] at %s\n',name,class(object),location); end
                            % call write for all members in the array
                            dumped_objects=0;
                            previous_msg='';
                            for n=1:numel(object)
                                dumped_objects=dumped_objects+h.write(object(n), [name '_' sprintf('%04d',n)],[location '/' name]);
                                if h.verbose
                                    msg = sprintf('Processed %d/%d', n, numel(object));
                                    fprintf([previous_msg, msg]);
                                    previous_msg = repmat(sprintf('\b'), 1, length(msg));
                                end
                            end
                            if h.verbose fprintf('\n'); end
                            
                            % group attributes
                            if dumped_objects
                                h5writeatt(h.filename,[location '/' name],'class',class(object));
                                h5writeatt(h.filename,[location '/' name],'name',name);
                                h5writeatt(h.filename,[location '/' name],'array',1);
                                h5writeatt(h.filename,[location '/' name],'size',size(object));
                            end
                        else
                            % here we process non-array UFF structures
                            switch class(object)
                                case {'uff.channel_data' 'uff.beamformed_data' 'uff.phantom'}
                                    if h.verbose fprintf('UFF: writting %s [%s] at %s\n',name,class(object),location); end
                            end
                            
                            % dump all fields in struct (or properties in class)
                            dumped_objects=0;
                            field_list = fieldnames(object);
                            eval(['mco = ?' class(object) ';']);
                            plist = mco.PropertyList;
                            for f=1:length(field_list)

                                
                                % check if the property is dependent
                                copy=false;
                                for p=1:length(plist)
                                    if strcmp(plist(p).Name,field_list{f})
                                        copy=~plist(p).Dependent;
                                        continue;
                                    end
                                end
                                
                                % if it isn't dependent or empty we write it
                                if copy
                                    prop=object.(field_list{f});
                                    if numel(prop)
                                        h.write(prop, field_list{f},[location '/' name]);
                                        dumped_objects=dumped_objects+1;
                                    end
                                end
                            end
                            
                            if dumped_objects>0
                                h5writeatt(h.filename,[location '/' name],'class',class(object));
                                h5writeatt(h.filename,[location '/' name],'name',name);
                                h5writeatt(h.filename,[location '/' name],'size',size(object));
                                h5writeatt(h.filename,[location '/' name],'array',0);
                            end
                        end
                    else warning(sprintf('Class %s not supported by UFF; skipping write.',class(object))); end
            end
        end
        
        
        function out=index(h,location,display)
            %% INDEX - Returns contains of a UFF location
            %
            %   This UFF method returns a list of the datasets and groups
            %   in the specified location. The location must be group. The
            %   display flag will plot the list on screen.
            %
            %   out = UFF.INDEX(location,display)
            %
            %   Example:
            %
            %       uff_file = uff('test.uff');
            %       list = uff_file.index('/',true);
            %
            %   See also UFF.READ, UFF.UFF, UFF.WRITE
            
            if nargin<2 location='/'; end
            if nargin<3 display=false; end
            if isempty(location) location='/'; end
            if isempty(display) display=false; end
            
            
            info=h5info(h.filename,location);
            if display fprintf('UFF: Contents of %s at %s\n',h.filename,location); end
            
            % groups
            groups=info.Groups;
            out={}; nn=1;
            for n=1:length(groups)
                out{nn}.location=groups(n).Name;
                out{nn}.name=[];
                out{nn}.class=[];
                out{nn}.size=[0 0];
                for m=1:length(groups(n).Attributes)
                    if strcmp(groups(n).Attributes(m).Name,'class') out{nn}.class=groups(n).Attributes(m).Value; end
                    if strcmp(groups(n).Attributes(m).Name,'name') out{nn}.name=groups(n).Attributes(m).Value; end
                    if strcmp(groups(n).Attributes(m).Name,'size') out{nn}.size=groups(n).Attributes(m).Value; end
                end
                if display fprintf('   - %s: %s [%s] size(%d,%d)\n',out{nn}.location, out{nn}.name, out{nn}.class, out{nn}.size); end
                nn=nn+1;
            end
            
            % datasets
            datasets=info.Datasets;
            for n=1:length(datasets)
                out{nn}.location=[location '/' datasets(n).Name];
                for m=1:length(datasets(n).Attributes)
                    if strcmp(datasets(n).Attributes(m).Name,'class') out{nn}.class=datasets(n).Attributes(m).Value; end
                    if strcmp(datasets(n).Attributes(m).Name,'name') out{nn}.name=datasets(n).Attributes(m).Value; end
                end
                if display fprintf('   - %s: %s [%s]\n',out{nn}.location, out{nn}.name, out{nn}.class); end
                nn=nn+1;
            end
            
        end
        
        
        function out = read(h, location)
            %% READ - Reads UFF dataset or group from location
            %
            %   This UFF method returns a variable containing the information
            %   in the specified location. The location can be a dataset or
            %   a group. Version v1.0.0.
            %
            %   out = UFF.READ(location)
            %
            %   Example:
            %
            %       uff_file = uff('test.uff');
            %       channel_data = uff_file('/channel_data');
            %
            %   See also UFF.WRITE, UFF.UFF, UFF.INDEX
            
            % checking version
            file_version=h5readatt(h.filename, '/','version');  % read file version
            file_version=file_version{1};                       % from cell to string
            file_version=file_version(int32(file_version)>0);   % removing 0's from 0-terminated strings
            
            switch file_version
                case h.version
                case 'v1.0.0'
                % here we will forward the read of obsolete versions to previous read functions
                otherwise
                    error(sprintf('Unsupported file version (%s). Current UFF version (%s). Please choose a new file instead.',file_version,h.version));
            end
            
            % check for root folder
            if nargin<2 || strcmp(location,'/')
                item=h.index('/');
                if length(item) out={}; end
                for n=1:length(item)
                    out{n}=h.read(item{n}.location);
                end
            else
                
                % checking name and class
                data_name=h5readatt(h.filename, location ,'name');
                class_name=h5readatt(h.filename, location ,'class');
                
                switch class_name
                    case {'double' 'single'}
                        if ~h5readatt(h.filename, location, 'complex')
                            out=h5read(h.filename, location);
                        else
                            out=h5read(h.filename, [ location '/real' ])+...
                                1i*h5read(h.filename, [ location '/imag' ]);
                        end
                    case 'char'
                        out=char(h5read(h.filename, location));
                    case 'uff.window'
                        out=uff.window(h5read(h.filename, location ));
                    case 'cell'
                        data_size=h5readatt(h.filename, location ,'size');
                        N=prod(data_size);
                        if(N>0)
                            item=h.index( location );
                            if length(item)~=N error('Size attribute does not match number of subgroups'); end
                            out={};
                            for n=1:N
                                out{n}=h.read(item{n}.location);
                            end
                            reshape(out,data_size.');
                        end
                    otherwise
                        % rest of UFF structures
                        if (findstr('uff.',class_name))
                            switch class_name
                                case {'uff.channel_data' 'uff.beamformed_data' 'uff.phantom'}
                                    if h.verbose fprintf('UFF: reading %s [%s]\n',data_name,class_name); end
%                                 case 'uff.wave'
%                                     if h.verbose fprintf('.'); end
                            end
                            data_size=h5readatt(h.filename, location ,'size');
                            N=prod(data_size);
                            if(N>1)
                                item=h.index( location );
                                if length(item)~=N error('Size attribute does not match number of subgroups'); end
                                if h.verbose fprintf('UFF: reading %s [%s]\n',data_name,class_name); end
                                previous_msg = '';
                                for n=1:N
                                    out(n)=h.read(item{n}.location);
                                    if h.verbose
                                        msg = sprintf('Processed %d/%d', n, N);
                                        fprintf([previous_msg, msg]);
                                        previous_msg = repmat(sprintf('\b'), 1, length(msg));
                                    end
                                end
                                if h.verbose fprintf('\n'); end
                                reshape(out,data_size.');
                            else
                                out=feval(class_name);
                                
                                % add properties
                                prop=h.index(location);
                                for m=1:length(prop)
                                    out.(prop{m}.name)=h.read(prop{m}.location);
                                end
                            end
                        else
                            warning(sprintf('Class %s not supported by UFF; skipping write.',class(value)));
                            out=[];
                        end
                end
            end
        end
    end
end