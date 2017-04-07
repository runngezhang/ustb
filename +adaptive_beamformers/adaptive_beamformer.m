classdef adaptive_beamformer
    %ADAPTIVE_BEAMFORMER defining the adaptive beamformers
    %
    %   All adaptive beamformeres have to be a subclass of this class
    %   See the coherence_factor as an example
    %
    %   authors: Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %
    %   $Date: 2017/04/07$
    
    %   All subclasses of the adaptive_beamformer will have access to these
    %   data parameters.
    properties (Access = protected)
        data_cube   % Delayed channel data from one sequence structured as
        % [z-samples,x-samples,channels]
        tx_apo      % Transmit apodization
        rx_apo      % Receive apodization
    end
    
    % These functions are used to set the protected variables
    methods (Access = public)
        function h = set_data_cube(h,data_cube)
            h.data_cube = data_cube;
        end
        function h = set_tx_apo(h,tx_apo)
            h.tx_apo = tx_apo;
        end
        function h = set_rx_apo(h,rx_apo)
            h.rx_apo = rx_apo;
        end
    end
    
    % All subclases of the adaptive_beamformer have to implement the "go"
    % method. Please see the coherence_factor for an example on how to
    % implement it
    methods (Abstract)
        go(obj)
    end
    
end

