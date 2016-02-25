function [sig] = ir_mex(h,r)
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
            sig=mex.cpwir(single(h.data), ...           % data
                single(r.scan.x.'),... 
                single(r.scan.z.'),...      % pixel positions (m)
                single([h.geom(:,1) h.geom(:,3)]),...   % probe geometry [x, z] (m,m)
                single(h.c0),...                % speed of sound (m/s)
                single(h.angle.'),...                   % angles of the plane waves (rad)
                single(h.transmit_apodization),...      % transmit aperture [pixels, channels]
                single(h.receive_apodization),...       % receive aperture [pixels, channels]
                single(h.sampling_frequency),...        % sampling frequency [Hz]
                single(h.initial_time),...              % initial time [s]
                single(0));
        case E.signal_format.IQ
            sig=mex.cpwir(single(h.data), ...       % data
                single(r.scan.x.'), ...
                single(r.scan.z.'),...      % pixel positions (m)
                single([h.geom(:,1) h.geom(:,3)]),...   % probe geometry [x, z] (m,m)
                single(h.c0),...                % speed of sound (m/s)
                single(h.angle).',...                   % angles of the plane waves (rad)
                single(h.transmit_apodization),...      % transmit aperture [pixels, channels]
                single(h.receive_apodization),...       % receive aperture [pixels, channels]
                single(h.sampling_frequency),...        % sampling frequency [Hz]
                single(h.initial_time),...              % initial time [s]
                single(h.modulation_frequency));        % modulation frequency [Hz]
        otherwise
            error('Unknown signal format!');
    end

end

