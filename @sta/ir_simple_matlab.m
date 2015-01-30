function sig = ir_simple_matlab(h,r)
%IR_SIMPLE_MATLAB    Image reconstruction implementation with simple matlab code
%
%   Syntax:
%   beamformed_signal = ir_simple_matlab(h,recons)
%       recons              RECONSTRUCTION class containing the specification of the reconstruction
%       beamformed_signal   Matrix containig the beamformed raw data 
%
%   See also RECONSTRUCTION

%   authors: Alfonso Rodriguez-Molares (alfonsom@ntnu.no)
%   $Date: 2015/01/28 $

    %% formats
    switch(h.format)
        case E.signal_format.RF
            w0=0;
        case E.signal_format.IQ
            w0=2*pi*h.modulation_frequency;
        otherwise
            error('Unknown signal format!');
    end
    
    %% beamforming
    sig=zeros(r.scan.pixels,h.frames);
    wb=waitbar(0,'STA Beamforming');
    for f=1:h.frames
        for ntx=1:h.channels
            waitbar(((f-1)*h.channels+ntx)/h.channels/h.frames);
            TF=sqrt((h.geom(ntx,1)-r.scan.x).^2+(h.geom(ntx,3)-r.scan.z).^2);
            for nrx=1:h.channels
                RF=sqrt((h.geom(nrx,1)-r.scan.x).^2+(h.geom(nrx,3)-r.scan.z).^2);
                delay=(RF+TF)/h.c0;
                phase_shift=exp(1i.*w0*delay); % <-- this is unnecessary in RF 
                sig(:,f)=sig(:,f)+phase_shift.*h.transmit_apodization(:,ntx).*h.receive_apodization(:,nrx).*interp1(h.time,h.data(:,nrx,ntx,f),delay,'linear',0);
            end
        end
    end
    close(wb);
    
end

