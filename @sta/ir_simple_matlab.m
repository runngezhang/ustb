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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RF format
        case E.signal_format.RF
            % beamforming
            sig=zeros(r.scan.pixels,h.F);
            wb=waitbar(0,'STA Beamforming, RF format');
            for f=1:h.F
                for nrx=1:h.M
                    waitbar(nrx/h.M);
                    RF=sqrt((h.geom(nrx,1)-r.scan.x).^2+(h.geom(nrx,3)-r.scan.z).^2);
                    for ntx=1:h.M
                        TF=sqrt((h.geom(ntx,1)-r.scan.x).^2+(h.geom(ntx,3)-r.scan.z).^2);
                        sig(:,f)=sig(:,f)+h.tx_apodization(:,ntx).*h.rx_apodization(:,nrx).*interp1(h.time,h.data(:,nrx,ntx,f),(RF+TF)/h.c0,'linear',0);
                    end
                end
            end
            close(wb);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IQ format
        case E.signal_format.IQ
            w0=2*pi*h.modulation_frequency;
            % beamforming
            sig=zeros(r.scan.pixels,h.F);
            wb=waitbar(0,'STA Beamforming, IQ format');
            for f=1:h.F
                for nrx=1:h.M
                    waitbar(nrx/h.M);
                    RF=sqrt((h.geom(nrx,1)-r.scan.x).^2+(h.geom(nrx,3)-r.scan.z).^2);
                    for ntx=1:h.M
                        TF=sqrt((h.geom(ntx,1)-r.scan.x).^2+(h.geom(ntx,3)-r.scan.z).^2);
                        delay=(RF+TF)/h.c0;
                        phase=w0*delay;
                        sig(:,f)=sig(:,f)+h.tx_apodization(:,ntx).*h.rx_apodization(:,nrx).*exp(1i.*phase).*interp1(h.time,h.data(:,nrx,ntx,f),delay,'linear',0);
                    end
                end
            end
            close(wb);
        otherwise
            error('Unknown signal format!');
    end
end

