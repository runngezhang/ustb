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
                
    % cheap workaround for matlab singular dimension problem
    unidimensional=false;
    if ndims(h.data)==2
        unidimensional=true;
        h.data=repmat(h.data,[1 1 1 2]);
    end
    
    switch(h.format)
        case E.signal_format.RF
            sig=mex.cpwir(h.data, ...           % data
                r.scan.x.', r.scan.z.',...      % pixel positions (m)
                [h.geom(:,1) h.geom(:,3)],...   % probe geometry [x, z] (m,m)
                double(h.c0),...                % speed of sound (m/s)
                h.angle.',...                   % angles of the plane waves (rad)
                h.transmit_apodization,...      % transmit aperture [pixels, channels]
                h.receive_apodization,...       % receive aperture [pixels, channels]
                h.sampling_frequency,...        % sampling frequency [Hz]
                h.initial_time,...              % initial time [s]
                0);
        case E.signal_format.IQ
            sig=mex.cpwir(double(h.data), ...       % data
                double(r.scan.x.'), double(r.scan.z.'),...      % pixel positions (m)
                double([h.geom(:,1) h.geom(:,3)]),...   % probe geometry [x, z] (m,m)
                double(h.c0),...                % speed of sound (m/s)
                double(h.angle).',...                   % angles of the plane waves (rad)
                double(h.transmit_apodization),...      % transmit aperture [pixels, channels]
                double(h.receive_apodization),...       % receive aperture [pixels, channels]
                double(h.sampling_frequency),...        % sampling frequency [Hz]
                double(h.initial_time),...              % initial time [s]
                double(h.modulation_frequency));        % modulation frequency [Hz]
        otherwise
            error('Unknown signal format!');
    end
    
    % cheap workaround for matlab singular dimension problem
    if unidimensional
        sig=sig(:,1);
    end
end

