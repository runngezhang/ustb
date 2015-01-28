function [sig] = ir_mex(h,r)
% function [sig] = ir_mex(h,r)

    disp('VS image reconstruction: mex implementation');

    % reshape pixel matrix specification
    rs_tx_apo=reshape(h.tx_apodization,[h.Nz*h.Nx h.N]); 
    rs_rx_apo=reshape(h.rx_apodization,[h.Nz*h.Nx h.M]);
                
    % Call mex
    switch(h.format)
        case E.signal_format.RF
            sig=mex.vsir_mex(h.data, ...       % data
                r.x(:).',r.z(:).',...           % pixel positions (m)
                [h.geom(:,1) h.geom(:,3)],...   % probe geometry [x, z] (m,m)
                h.c0,...                        % speed of sound (m/s)
                [h.source(:,1) h.source(:,3)],...% virtual sources location [x, z] (m,m)
                rs_tx_apo,...                   % transmit aperture [pixels, channels]
                rs_rx_apo,...                   % receive aperture [pixels, channels]
                h.Fs,...                        % sampling frequency [Hz]
                h.time(1));                     % initial time [s]
        case E.signal_format.IQ
            sig=mex.vsir_mex(h.data, ...       % data
                r.x(:).',r.z(:).',...           % pixel positions (m)
                [h.geom(:,1) h.geom(:,3)],...   % probe geometry [x, z] (m,m)
                h.c0,...                        % speed of sound (m/s)
                [h.source(:,1) h.source(:,3)],...% virtual sources location [x, z] (m,m)
                rs_tx_apo,...                   % transmit aperture [pixels, channels]
                rs_rx_apo,...                   % receive aperture [pixels, channels]
                h.Fs,...                        % sampling frequency [Hz]
                h.time(1),...                   % initial time [s]
                h.modulation_frequency);        % modulation frequency [Hz]
        otherwise
            error('Unknown signal format!');
    end

    % reshape beamformed image
    sig=reshape(sig,[h.Nz h.Nx h.F]);
    
end

