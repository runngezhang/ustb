function [sig] = ir_simple_matlab(h,r)
% function [sig] = ir_simple_matlab(h)

    disp('CPW image reconstruction: simple matlab implementation');
    disp('alfonso.r.molares@ntnu.no, 23/01/2015 ');

    %% formats
    switch(h.format)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RF format
        case E.signal_format.RF
            
            % beamforming
            sig=zeros(h.Nz,h.Nx,h.F);
            wb=waitbar(0,'CPW Beamforming, RF format');
            for f=1:h.F
                for na=1:h.N
                    waitbar(na/h.N);
                    Da=r.z*cos(h.angle(na))+r.x*sin(h.angle(na));
                    for nrx=1:h.M
                        RF=sqrt((h.geom(nrx,1)-r.x).^2+(h.geom(nrx,3)-r.z).^2);
                        sig(:,:,f)=sig(:,:,f)+h.tx_apodization(:,:,na).*h.rx_apodization(:,:,nrx).*interp1(h.time,h.data(:,nrx,na,f),(Da+RF)/h.c0,'linear',0);
                    end
                end
            end
            close(wb);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IQ format
        case E.signal_format.IQ
            w0=2*pi*h.modulation_frequency;
            
            % beamforming
            sig=zeros(h.Nz,h.Nx,h.F);
            wb=waitbar(0,'CPW Beamforming, IQ format');
            for f=1:h.F
                for na=1:h.N
                    waitbar(na/h.N);
                    Da=r.z*cos(h.angle(na))+r.x*sin(h.angle(na));
                    for nrx=1:h.M
                        RF=sqrt((h.geom(nrx,1)-r.x).^2+(h.geom(nrx,3)-r.z).^2);
                        delay=(Da+RF)/h.c0;
                        phase=w0*delay;
                        sig(:,:,f)=sig(:,:,f)+h.tx_apodization(:,:,na).*h.rx_apodization(:,:,nrx).*exp(1i.*phase).*interp1(h.time,h.data(:,nrx,na,f),delay,'linear',0);
                    end
                end
            end
            close(wb);
            
        otherwise
            error('Unknown signal format!');
    end
end

