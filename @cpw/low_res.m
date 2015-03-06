function [sig] = low_res(h,r)
%LOW_RES    Provide low resolution images
%
%   Syntax:
%   beamformed_signal = low_res(h,recons)
%       recons              RECONSTRUCTION class containing the specification of the reconstruction
%       beamformed_signal   Matrix containig the beamformed raw data 
%
%   See also RECONSTRUCTION

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/03/06 $


    % initial time -> performance test
    tic; 

    % loop over orientations
    total_data=zeros(size(r.scan.x_matrix,1),size(r.scan.x_matrix,2),length(r.orientation),size(h.data,4));
    assert(length(r.orientation)==1,'Calculation of low resolution images only works for single orientation. To be clarified.')

    % precompute transmit and receive apodization
    xT=r.scan.x*ones(1,h.firings)-r.scan.z*tan(h.angle.'); % position of equivalent receive element -> Alfonso's equation 
    h.transmit_apodization = r.calculate_apodization(r.orientation.transmit_beam,xT);
    xR=ones(r.scan.pixels,1)*(h.geom(:,1).');                  % position of receive element
    h.receive_apodization = r.calculate_apodization(r.orientation.receive_beam,xR);

    switch(h.format)
        case E.signal_format.RF
            sig=mex.cpwlr(h.data, ...       % data
                r.scan.x.', r.scan.z.',...      % pixel positions (m)
                [h.geom(:,1) h.geom(:,3)],...   % probe geometry [x, z] (m,m)
                h.c0,...                        % speed of sound (m/s)
                h.angle.',...                   % angles of the plane waves (rad)
                h.transmit_apodization,...      % transmit aperture [pixels, channels]
                h.receive_apodization,...       % receive aperture [pixels, channels]
                h.sampling_frequency,...        % sampling frequency [Hz]
                h.initial_time);                % initial time [s]
        case E.signal_format.IQ
            sig=mex.cpwlr(h.data, ...       % data
                r.scan.x.', r.scan.z.',...      % pixel positions (m)
                [h.geom(:,1) h.geom(:,3)],...   % probe geometry [x, z] (m,m)
                h.c0,...                        % speed of sound (m/s)
                h.angle.',...                   % angles of the plane waves (rad)
                h.transmit_apodization,...      % transmit aperture [pixels, channels]
                h.receive_apodization,...       % receive aperture [pixels, channels]
                h.sampling_frequency,...        % sampling frequency [Hz]
                h.initial_time,...              % initial time [s]
                h.modulation_frequency);        % modulation frequency [Hz]
        otherwise
            error('Unknown signal format!');
    end

end

