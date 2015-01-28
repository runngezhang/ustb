classdef dataset < handle
    % dataset Superclass for different kinds of datasets
    %   dataset class contains all the information and methods that are
    %   common to all kinds of datasets. Datasets could be Synthetic transmit
    %   aperture (sta), Coherent plane wave (cpw), Virtual source (vs), 
    %   and even image reconstructions (ir).
    
    properties  (SetAccess = public)
        name            % String containing the name of the dataset
        creation_date   % String containing the date the dataset class was created
        format = E.signal_format.IQ  % Signal_format RF or IQ (enumerations.signal_format.RF, default=enumerations.signal_format.IQ)
        geom            % matrix M x 3 containing probe geometry [x, y, z] (m)
        data            % Matrix containing the numerical data. For acquisition
                        % datasets the matrix dimensions are [time_samples, channels, firings, frames]
                        % For image reconstruction datasets the matrix
                        % dimensions are [pixels, frames]
        time            % vector containing fast time (s)
        c0              % value conatining  the reference speed of sound (m/s)
        modulation_frequency % value conatining the modulation frequency (Hz), only required for IQ format
    end
    properties  (SetAccess = protected)   
        F               % number of frames in the dataset        
        M               % number of elements in the transducer
        Nx              % Number of pixel in the x direction (to be removed)
        Nz              % Number of pixel in the z direction (to be removed)
        t0              % initial time (s)
        Fs              % sampling frequency (Hz)
        tx_apodization  % matrix containing the apodization used for transmit
        rx_apodization  % matrix containing the apodization used for transmit
        elapsed_time    % elapsed time, for performance control
    end
    
    methods  % formating methods
        function set.time(h,input_time)
            if(size(input_time,1)<size(input_time,2))
                error('The time vector must be a column vector!');
            else
                h.time=input_time;
            end
        end
    end
    
    methods (Access = public)
        function h = dataset(name)
            % The constructor only requires the name of the dataset. It
            % automatically sets the creation date.
            h.name=name;
            t_now=now;
            h.creation_date=[date sprintf('-%d-%d-%d',hour(t_now),minute(t_now),round(second(t_now)))];
        end
        
        function [apo]= linear_apodization(h,beam,r)
            if beam.apodization==E.apodization_type.none
                apo=ones(h.Nz,h.Nx,h.M);
            else
                apo=zeros(h.Nz,h.Nx,h.M);
                Aperture=r.z./beam.FN;
                for n=1:h.M
                    xd=abs(h.geom(n,1)-r.x+r.z*tan(beam.steer_angle));
                    switch(beam.apodization)
                        case E.apodization_type.boxcar
                            apo(:,:,n)=double(xd<Aperture/2); 
                        case E.apodization_type.hanning
                            apo(:,:,n)=double(xd<Aperture/2).*(0.5 + 0.5*cos(2*pi*xd./Aperture)); 
                        case E.apodization_type.tukey25
                            roll=0.25;
                            apo(:,:,n)=(xd<(Aperture/2*(1-roll))) + (xd>(Aperture/2*(1-roll))).*(xd<(Aperture/2)).*0.5.*(1+cos(2*pi/roll*(xd/Aperture-roll/2-1/2)));                               
                        case E.apodization_type.tukey50
                            roll=0.5;
                            apo(:,:,n)=(xd<(Aperture/2*(1-roll))) + (xd>(Aperture/2*(1-roll))).*(xd<(Aperture/2)).*0.5.*(1+cos(2*pi/roll*(xd/Aperture-roll/2-1/2)));                               
                        case E.apodization_type.tukey80
                            roll=0.8;
                            apo(:,:,n)=(xd<(Aperture/2*(1-roll))) + (xd>(Aperture/2*(1-roll))).*(xd<(Aperture/2)).*0.5.*(1+cos(2*pi/roll*(xd/Aperture-roll/2-1/2)));                               
                        otherwise
                            error('Unknown apodization type!');
                    end
                end
            end
            % implementation of Bastian's damping
            if(beam.damping)
                n0=1:(beam.damping+1);                                  
                mask=(n0/(beam.damping+1)).^(1./beam.damping_order);
                for nn=n0(1):n0(end) % have to think on how to vectorise this
                    apo(:,:,nn)=apo(:,:,nn).*mask(nn);
                    apo(:,:,end-nn+1)=apo(:,:,end-nn+1).*mask(nn);
                end
            end
        end
        
        function [im]=envelope(h,sig)
            % function [im]=envelope(h,sig)
            switch(h.format)
                case E.signal_format.RF
                    warning('Computing Hilbert transform on pixel data. Envelope may be inaccurate for low mesh resolution.');
                    im=zeros(h.Nz,h.Nx,h.F);
                    for f=1:h.F
                        for nx=1:h.Nx
                            im(:,nx,f)=abs(hilbert(sig(:,nx,f)));
                        end
                    end
                case E.signal_format.IQ
                    im=abs(sig);
                otherwise
                    error('Unknown signal format!');
            end
        end
        
        function show(h,r,input_im,dynamic_range)
            % Ploting image reconstruction
            im=20*log10(input_im./max(input_im(:)));
            figure; set(gca,'fontsize',16);
            for f=1:h.F
                pcolor(r.x*1e3,r.z*1e3,im(:,:,f)); shading flat; axis equal; axis tight; colormap gray; caxis([-dynamic_range 0]);colorbar;
                xlabel('x [mm]');
                ylabel('z [mm]');
                set(gca,'YDir','reverse');
                set(gca,'fontsize',16);
                axis([min(r.x(:)) max(r.x(:)) min(r.z(:)) max(r.z(:))]*1e3);
                title(sprintf('%s (%s, %0.2fs/f)',char(h.name),char(h.format),h.elapsed_time/h.F)); 
                set(gcf,'Position',[680   618   640   480]);
                drawnow;
                pause(0.05);
            end
        end
        
        function [sig]=launch_implementation(h,r,imp)
            % Launches appropriate implementation
            switch(imp)
                case E.implementation.simple_matlab
                    sig=h.ir_simple_matlab(r);
                case E.implementation.optimized_matlab
                    sig=h.ir_optimized_matlab(r);
                case E.implementation.mex
                    sig=h.ir_mex(r);
                case E.implementation.mex_gpu
                    sig=h.ir_mex_gpu(r);
                otherwise
                    error('Selected implementation is not supported');
            end
        end
    end
    
    methods 
        function sig=ir_simple_matlab(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Simple matlab: The implementation is not available yet!');
            sig=[];
        end
        function sig=ir_optimized_matlab(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Optimazed matlab: The implementation is not available yet!');
            sig=[];
        end
        function sig=ir_mex(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Mex: The implementation is not available yet!');
            sig=[];
        end
        function sig=ir_mex_gpu(h,r)
            % Simple Matlab implementation to be overloaded by subclasses
            error('Mex gpu: The implementation is not available yet!');
            sig=[];
        end
    end
end

