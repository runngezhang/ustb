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

    switch(h.format)
        case E.signal_format.RF
            sig=mex.stair_mex(h.data, ...       % data
                r.scan.x.', r.scan.z.',...      % pixel positions (m)
                [h.geom(:,1) h.geom(:,3)],...   % probe geometry [x, z] (m,m)
                h.c0,...                        % speed of sound (m/s)
                h.transmit_apodization,...      % transmit aperture [pixels, channels]
                h.receive_apodization,...       % receive aperture [pixels, channels]
                h.sampling_frequency,...        % sampling frequency [Hz]
                h.initial_time);                % initial time [s]
        case E.signal_format.IQ
            sig=mex.stair_mex(h.data, ...       % data
                r.scan.x.', r.scan.z.',...      % pixel positions (m)
                [h.geom(:,1) h.geom(:,3)],...   % probe geometry [x, z] (m,m)
                h.c0,...                        % speed of sound (m/s)
                h.transmit_apodization,...      % transmit aperture [pixels, channels]
                h.receive_apodization,...       % receive aperture [pixels, channels]
                h.sampling_frequency,...        % sampling frequency [Hz]
                h.initial_time,...              % initial time [s]
                h.modulation_frequency);        % modulation frequency [Hz]
        otherwise
            error('Unknown signal format!');
    end
    
end

