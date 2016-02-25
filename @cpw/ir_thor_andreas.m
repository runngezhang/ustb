function [sig] = ir_thor_andreas(h,r)
%IR_MEX    Thor-Andreas beamforming code, mex compilation
%
%   Syntax:
%   beamformed_signal = ir_mex(h,recons)
%       recons              RECONSTRUCTION class containing the specification of the reconstruction
%       beamformed_signal   Matrix containig the beamformed raw data 
%
%   See also RECONSTRUCTION

%   authors: Thor-Andreas Tangen ??

    % Parameters
    p           = {};
    p.c         = single(h.c0);
    p.pitch     = single(h.geom(2,1)-h.geom(1,1));
    p.dx        = p.pitch;
    p.t0        = single(h.time(1));
    p.kerf      = single(0); % needed?
    p.tx_angle  = single(h.angle);
    p.rx_angle  = single(zeros(length(h.angle),1));
    p.fs_in     = single(h.sampling_frequency);
    p.fs_out    = single(p.fs_in);
    p.f_demod   = single(h.modulation_frequency);
    p.FN        = single(1); % needed?

%     keyboard
%     % coincident z-axis condition
%     ta_z_axis=h.time*h.c0/2;
%     n_end=2*r.scan.z_axis(end)/h.c0*h.sampling_frequency;
%     ta_z_axis=ta_z_axis(1:n_end);

    % call mex
    sig=zeros(r.scan.pixels,h.frames);
    for n=1:h.frames
        sig_temporal = mex.planewave_beamforming2(single(h.data),p,single(r.scan.x_axis),single(r.scan.z_axis),single(reshape(h.receive_apodization,[r.scan.Nz r.scan.Nx h.channels])));  % beamforming procedure
        sig_temporal(isnan(sig_temporal)) = 0;
        
        % compound
        sig_temporal = sum(sig_temporal,3);

        % undo the distance correction factor
        ta_z_axis=h.time*h.c0/2;
        corr_factor=exp(j*2*pi*h.modulation_frequency*h.time);
        sig_temporal=bsxfun(@times,sig_temporal,corr_factor(1:length(r.scan.z_axis)));
        
        % reshape
        sig(:,n) = sig_temporal(:);
    end
    
    % correct TA distance factor
    %sig = sig.*
end

