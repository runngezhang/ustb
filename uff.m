classdef uff < handle
    %UFF   Ultrasound File Format (UFF) superclass
    %
    %   See also CHANNEL_DATA, BEAMFORMED_DATA
    
    %   authors: Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
    %
    %   $Last updated: 2017/06/07$
    
    %% Logistics parameters
    properties  (SetAccess = public)
        name={}                         % name of the dataset
        reference={}                    % reference to the publication where it was used/acquired
        author={}                       % contact of the authors
        version={}                      % version of the dataset
        info={}                         % other information
    end
    
    %% copy
    methods (Access = public)
        function copy(h,object)
            %COPY    Copy the values from another channel_data
            %
            %   Syntax:
            %   COPY(object)
            %       object       Instance of a channel_data class
            %
            %   See also SCAN, WAVE, SOURCE
            
            assert(isa(object,class(h)),'Class of the input object is not identical');
            
            % we copy all non-dependent & public properties
            list_properties=properties(object);
            for n=1:numel(list_properties)
                property_name=list_properties{n};
                mp = findprop(h,property_name);
                if strcmp(mp.GetAccess,'public')&&~mp.Dependent
                    eval(sprintf('h.%s = object.%s;',property_name,property_name));
                end
            end
        end
    end
    
    methods
        
        function  out  = objname(h)
            %% OBJNAME Gives object name
            out = evalin('caller','inputname(1)');
        end
        
        function out = hash(h)
            %% HASH Gives hash of the whole uff class
            
            % loop over all non-dependent & public properties
            str=[];
            list_properties=properties(h);
            for n=1:numel(list_properties)
                property_name=list_properties{n};
                mp = findprop(h,property_name);
                if strcmp(mp.GetAccess,'public')&&~mp.Dependent
                    %fprintf(1,'%s -> %s\n',property_name,tools.hash(h.(property_name)));
                    str=[str;tools.hash(h.(property_name))];
                end
            end
            
            out=tools.hash(str);
        end
    end
    
    %% Public UFF file write/read
    methods (Access = public)
        function write(h, filename, name, location, verbose)
            %% WRITE  Writes object into location
            %
            %   This UFF method writes the object into the specified location
            %   with the provided name.
            %
            %   UFF_OBJECT.WRITE(filename,name,location,verbose)
            %
            %   Parameters:
            %       filename    Name and path to the UFF file
            %       name        Name of the object within the file
            %       location    Location within the UFF file
            %       verbose     Flag to get text messages
            %
            %   Example:
            %       channel_data = uff.channel_data();
            %       channel_data.write('test.uff');
            %
            %   See also UFF.READ, UFF.INDEX
            
            if nargin<3||isempty(name) name=h.objname; end
            if nargin<4 location=[]; end
            if nargin<5||isempty(verbose) verbose=true; end
            
            % we are good to write the object
            uff.write_object(filename, h, name, location, verbose);            
        end
        
         function read(h, filename, location, verbose)
            %% READ  Reads object from location
            %
            %   This UFF method read the UFF object from the specified location
            %   with the provided name. 
            %
            %   UFF_OBJECT.READ(filename, location, verbose)
            %
            %   Parameters:
            %       filename    Name and path to the UFF file
            %       location    Location within the UFF file
            %       verbose     Flag to get text messages
            %
            %   Example:
            %       channel_data = uff.channel_data();
            %       channel_data.read('test.uff','/channel_data');
            %
            %   See also UFF.WRITE, UFF.INDEX
            
            if nargin<4 verbose=true; end
            
            object=uff.read_object(filename, location, verbose);
            if numel(object)==1 h.copy(object); 
            else
                error('UFF: Trying to copy an array of objects into an UFF class. Use intead: a = uff().read(filename,location)');
            end
        end
    end
    
    %% Private UFF file read/write
    methods (Access = private)
      
        
        
       
        
    end
    
    %% display methods
    methods
        function print_authorship(h)
            out_name = textwrap([],h.name,50);
            fprintf('Name: \t\t %s \n',out_name{1});
            for i = 2:numel(out_name)
                fprintf('\t\t %s \n',out_name{i});
            end
            fprintf('Reference: \t %s \n',h.reference{:});
            fprintf('Author(s): ');
            for i = 1:numel(h.author)
                if i == 1
                    fprintf('\t %s \n',h.author{i});
                else
                    fprintf('\t\t %s \n',h.author{i});
                end
            end
            fprintf('Version: \t %s \n',h.version{:});
        end
    end
end