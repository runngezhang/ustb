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
        version = 'v1.0.0';
    end
    
    %% constructor
    methods (Access = public)
        function h=uff(filename,mode)
            
            if nargin<1
                error('UFF: Missing filename');
            end
            if nargin<2
                mode='read';
            end
            
            switch mode
                case 'write'
                    attr_details.Name = 'version';
                    attr_details.AttachedTo = '/';
                    attr_details.AttachType = 'group';
                    hdf5write(filename, attr_details, h.version);
                case {'read', 'append'}
                otherwise
                    error(sprintf('Mode %s not supported',mode));
            end
            
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
%         function dumped_props=write(h,location,obj,name)
%             %% dumped_props=write(h,location,obj,name)
%             %
%             % Writes the object/structure obj at location with name
%             if(nargin<3)
%                 error('Missing object to write on UFF file');
%             end
%             if(nargin<4)
%                 error('Missing object name');
%             end
%             
%             % dump all fields in struct
%             dumped_props=0;
%             field_list = fieldnames(obj);
%             eval(['mco = ?' class(obj) ';']);
%             plist = mco.PropertyList;
%             for f=1:length(field_list)
%                 prop=obj.(field_list{f});
%                 
%                 copy=false;
%                 for p=1:length(plist)
%                     if strcmp(plist(p).Name,field_list{f})
%                         copy=~plist(p).Dependent;
%                         continue;
%                     end
%                 end
%                 
%                 if copy
%                     if numel(prop)
%                         h.write_field([location '/' name], field_list{f}, prop);
%                         dumped_props=dumped_props+1;
%                     end
%                 end
%             end
%             
%             if dumped_props>0
%                 h5writeatt(h.filename,[location '/' name],'class',class(obj));
%                 h5writeatt(h.filename,[location '/' name],'name',name);
%                 h5writeatt(h.filename,[location '/' name],'size',size(obj));
%             end
%         end
%         
%         function write_field(h, location, name, value)
%             %% out=write_field(h, location, name, value)
%             %
%             % Writes the field
%             switch class(value)
%                 case {'double' 'single'}
%                     if isreal(value)
%                         h5create(h.filename,[location '/' name], size(value), 'Datatype', 'single', 'ChunkSize',size(value));
%                         h5write(h.filename,[location '/' name], single(value));
%                         h5writeatt(h.filename,[location '/' name],'class',class(value));
%                         h5writeatt(h.filename,[location '/' name],'name',name);
%                         h5writeatt(h.filename,[location '/' name],'imaginary',0);
%                         h5writeatt(h.filename,[location '/' name],'complex',0);
%                     else
%                         % real
%                         h5create(h.filename,[location '/' name '/real'], size(value), 'Datatype', 'single', 'ChunkSize',size(value));
%                         h5write(h.filename,[location '/' name '/real'], single(real(value)));
%                         h5writeatt(h.filename,[location '/' name '/real'],'class',class(value));
%                         h5writeatt(h.filename,[location '/' name '/real'],'name',name);
%                         h5writeatt(h.filename,[location '/' name '/real'],'imaginary',0);
%                         
%                         % imag
%                         h5create(h.filename,[location '/' name '/imag'], size(value), 'Datatype', 'single', 'ChunkSize',size(value));
%                         h5write(h.filename,[location '/' name '/imag'], single(imag(value)));
%                         h5writeatt(h.filename,[location '/' name '/imag'],'class',class(value));
%                         h5writeatt(h.filename,[location '/' name '/imag'],'name',name);
%                         h5writeatt(h.filename,[location '/' name '/imag'],'imaginary',1);
%                         
%                         % group attributes
%                         h5writeatt(h.filename,[location '/' name],'class',class(value));
%                         h5writeatt(h.filename,[location '/' name],'name',name);
%                         h5writeatt(h.filename,[location '/' name],'complex',1);
%                     end
%                 case {'char' 'uff.window'}
%                     h5create(h.filename,[location '/' name], size(value), 'Datatype', 'single', 'ChunkSize',size(value));
%                     h5write(h.filename,[location '/' name], single(value));
%                     h5writeatt(h.filename,[location '/' name],'class',class(value));
%                     h5writeatt(h.filename,[location '/' name],'name',name);
%                 otherwise
%                     % UFF structures
%                     if (findstr('uff.',class(value)))
%                         if numel(value)>1
%                             dumped_props=0;
%                             for n=1:numel(value)
%                                 dumped_props=dumped_props+h.write([location '/' name], value(n), [name '_' sprintf('%04d',n)]);
%                             end
%                             
%                             % group attributes
%                             if dumped_props
%                                 h5writeatt(h.filename,[location '/' name],'class',class(value));
%                                 h5writeatt(h.filename,[location '/' name],'name',name);
%                                 h5writeatt(h.filename,[location '/' name],'array',1);
%                                 h5writeatt(h.filename,[location '/' name],'size',size(value));
%                             end
%                         else
%                             dumped_props=h.write(location, value, name);
%                             if dumped_props
%                                 h5writeatt(h.filename,[location '/' name],'array',0);
%                                 h5writeatt(h.filename,[location '/' name],'size',size(value));
%                             end
%                         end
%                     else warning(sprintf('Class %s not supported by UFF; skipping write.',class(value))); end
%             end
%         end
        
        function dumped_objects=write(h,location,object,name)
            %% WRITE- Writes object into location
            % 
            %   This UFF method writes the object into the specified location 
            %   with the provided name. A UFF object must be defined.
            %
            %   dumped_objects=UFF.WRITE(location,object,name)
            %
            %   Example:
            %
            %       channel_data = uff.channel_data();
            %       uff_file = uff('test.uff');
            %       dumped_objects = uff_file.write('/',channel_data,'chn01');
            %
            %   See also UFF.READ, UFF.UFF, UFF.INDEX

            if(nargin<3)||isempty(object)
                error('Missing object to write on UFF file');
            end
            if(nargin<4)||isempty(name)
                error('Missing object name');
            end
            
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
                case {'char' 'uff.window'}
                    h5create(h.filename,[location '/' name], size(object), 'Datatype', 'single', 'ChunkSize',size(object));
                    h5write(h.filename,[location '/' name], single(object));
                    h5writeatt(h.filename,[location '/' name],'class',class(object));
                    h5writeatt(h.filename,[location '/' name],'name',name);
                    dumped_objects=1;
                otherwise
                    % UFF structures
                    if (findstr('uff.',class(object)))
                        if numel(object)>1
                            
                            % call write for all members in the array
                            dumped_objects=0;
                            for n=1:numel(object)
                                dumped_objects=dumped_objects+h.write([location '/' name], object(n), [name '_' sprintf('%04d',n)]);
                            end
                            
                            % group attributes
                            if dumped_objects
                                h5writeatt(h.filename,[location '/' name],'class',class(object));
                                h5writeatt(h.filename,[location '/' name],'name',name);
                                h5writeatt(h.filename,[location '/' name],'array',1);
                                h5writeatt(h.filename,[location '/' name],'size',size(object));
                            end
                        else
                            % here we process the single UFF structures
                            %dumped_objects=h.write(location, object, name);
               
                            % dump all fields in struct (or properties in class)
                            dumped_objects=0;
                            field_list = fieldnames(object);
                            eval(['mco = ?' class(object) ';']);
                            plist = mco.PropertyList;
                            for f=1:length(field_list)
                                prop=object.(field_list{f});

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
                                    if numel(prop)
                                        h.write([location '/' name], prop, field_list{f});
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
            if display fprintf('Contents of %s at %s\n',h.filename,location); end
            
            % groups
            groups=info.Groups;
            out={}; nn=1;
            for n=1:length(groups)
                out{nn}.location=groups(n).Name;
                for m=1:length(groups(n).Attributes)
                    if strcmp(groups(n).Attributes(m).Name,'class') out{nn}.class=groups(n).Attributes(m).Value; end
                    if strcmp(groups(n).Attributes(m).Name,'name') out{nn}.name=groups(n).Attributes(m).Value; end
                end
                if display fprintf(' -- %s: %s [%s]\n',out{nn}.location, out{nn}.name, out{nn}.class); end
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
                if display fprintf(' -- %s: %s [%s]\n',out{nn}.location, out{nn}.name, out{nn}.class); end
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
            version='v1.0.0';                                   % this code version
            file_version=h5readatt(h.filename, '/','version');  % read file version
            file_version=file_version{1};                       % from cell to string
            file_version=file_version(int32(file_version)>0);   % removing 0's from 0-terminated strings
            assert(strcmp(version,file_version),sprintf('The file version (%s) and the READ function version (%s) do not match.',file_version,version));
            
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
                    out=h5read(h.filename, location)
                case 'uff.window'
                    out=uff.window(h5read(h.filename, location ));
                otherwise
                    % rest of UFF structures
                    if (findstr('uff.',class_name))
                        data_size=h5readatt(h.filename, location ,'size');
                        N=prod(data_size);
                        if(N>1)
                            item=h.index( location );
                            if length(item)~=N error('Size attribute does not match number of subgroups'); end
                            for n=1:N
                                out(n)=h.read(item{n}.location);
                            end
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