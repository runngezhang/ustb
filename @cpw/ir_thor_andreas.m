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
    p.fs_in     = single(h.sampling_frequency);
    p.fs_out    = single(p.fs_in);
    p.f_demod   = single(h.modulation_frequency);
    p.FN        = single(1); % needed?
    p.tx_angle  = single(repmat(h.angle,[h.frames]));
    p.rx_angle  = single(zeros(h.firings*h.frames,1));
    
    % call mex
    sig=zeros(r.scan.pixels,h.firings,h.frames);

    sig_temporal = mex.planewave_beamforming2(single(h.data(:,:,:)),p,single(r.scan.x_axis),single(r.scan.z_axis),single(reshape(h.receive_apodization,[r.scan.Nz r.scan.Nx h.channels])));  % beamforming procedure
    sig_temporal(isnan(sig_temporal)) = 0;

    % undo the distance correction factor
    ta_z_axis=h.time*h.c0/2;
    corr_factor=exp(j*2*pi*h.modulation_frequency*h.time);
    sig_temporal=bsxfun(@times,sig_temporal,corr_factor(1:r.scan.Nz));
        
    % reshape
    sig = reshape(sig_temporal,[r.scan.pixels h.firings h.frames]);
   

end

