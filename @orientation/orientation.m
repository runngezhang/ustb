classdef orientation < handle
%ORIENTATION   View class
%
%   See also RECONSTRUCTION, BEAM, SCAN
    properties  (SetAccess = public)
        transmit_beam=beam()        % BEAM class defining transmit beam
        receive_beam=beam()         % BEAM class defining receive beam            
    end
    
    %% constructor
    methods (Access = public)
        function h = orientation(tx_beam,rx_beam)
            %ORIENTATION    Constructor of the ORIENTATION class.
            %
            %   Syntax:
            %   ORIENTATION(transmit_beam, receive_beam) 
            %       transmit_beam           Definition of the transmit beam
            %       receive_beam            Definition of the receive beam
            %
            %   See also BEAM, LINEAR_SCAN
                      
            if exist('tx_beam') 
                h.transmit_beam=tx_beam;
            end

            if exist('rx_beam') 
                h.receive_beam=rx_beam;
            end
        end
    end
    
    %% copy
    methods (Access = public)
        function copy(h,object)
            %COPY    Copy the values from another ORIENTATION object
            %
            %   Syntax:
            %   COPY(object) 
            %       object       Instance of a ORIENTATION class
            %
            %   See also BEAM, RECONSTRUCTION
            assert(isa(object,class(h)),'Class of the input object is not identical'); 

            h.transmit_beam=object.transmit_beam;   % beam is not a handle class thus we can copy it
            h.receive_beam=object.receive_beam;     % beam is not a handle class thus we can copy it       
        end
    end
    
    %% huff 
    methods (Access = public)
        function huff_read(h,filename,location)
            %HUFF_READ    Read all the information of the orientation to a group in a HDF5
            %
            %   Syntax:
            %   HUFF_READ(file_name,location)
            %       file_name                  Name of the hdf5 file
            %       location                   Name of the group destination
            %
            %   See also ORIENTATION

            h.transmit_beam=h.transmit_beam.huff_read(filename,[location '\transmit_beam']);
            h.receive_beam=h.receive_beam.huff_read(filename,[location '\receive_beam']);
        end
        
         function huff_write(h,filename,location)
            %HUFF_WRITE    Dumps all the information of the orientation to a group in a HDF5
            %
            %   Syntax:
            %   HUFF_WRITE(file_name,location)
            %       file_name                  Name of the hdf5 file
            %       location                   Name of the group destination
            %
            %   See also ORIENTATION
         
            h.transmit_beam.huff_write(filename,[location '/transmit_beam']);
            h.receive_beam.huff_write(filename,[location '/receive_beam']);
         end
    end
    
    
end

