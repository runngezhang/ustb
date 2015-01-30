function sig = ir_mex(h,r)
%IR_MEX    Image reconstruction implementation with mex file
%
%   Syntax:
%   beamformed_signal = ir_mex(h,recons)
%       recons              RECONSTRUCTION class containing the specification of the reconstruction
%       beamformed_signal   Matrix containig the beamformed raw data 
%
%   See also RECONSTRUCTION

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/01/28 $
                
    % Call mex
    switch(h.format)
        case E.signal_format.RF
            sig=mex.vsir_mex(h.data, ...        % data
                r.scan.x.', r.scan.z.',...      % pixel positions (m)
                [h.geom(:,1) h.geom(:,3)],...   % probe geometry [x, z] (m,m)
                h.c0,...                        % speed of sound (m/s)
                [h.source(:,1) h.source(:,3)],...% virtual sources location [x, z] (m,m)
                h.tx_apodization,...            % transmit aperture [pixels, channels]
                h.rx_apodization,...            % receive aperture [pixels, channels]
                h.Fs,...                        % sampling frequency [Hz]
                h.time(1));                     % initial time [s]
        case E.signal_format.IQ
            sig=mex.vsir_mex(h.data, ...        % data
                r.scan.x.', r.scan.z.',...      % pixel positions (m)
                [h.geom(:,1) h.geom(:,3)],...   % probe geometry [x, z] (m,m)
                h.c0,...                        % speed of sound (m/s)
                [h.source(:,1) h.source(:,3)],...% virtual sources location [x, z] (m,m)
                h.tx_apodization,...            % transmit aperture [pixels, channels]
                h.rx_apodization,...            % receive aperture [pixels, channels]
                h.Fs,...                        % sampling frequency [Hz]
                h.time(1),...                   % initial time [s]
                h.modulation_frequency);        % modulation frequency [Hz]
        otherwise
            error('Unknown signal format!');
    end
    
end

