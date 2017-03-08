classdef verasonics < handle
    %ALPINION   Class for reading Verasonics research scanner data to USTB
    
    %   authors:    Alfonso Rodriques-Morales
    %               Ole Marius Hoel Rindal (olemarius@olemarius.net)
    %   $Date: 2017/02/17$
    
    properties (SetAccess = public)
        % Verasonics objects and structs
        Trans                  % Verasonics Transducer object
        TW                     % Verasonics Transmit Waveform object
        RcvData                % Verasonics Receive Data buffer, the channeldata
        Resource               % Verasonics Resource object, a few global parameters
        Receive                % Verasonics Receive object, defining receive parameters
        TX                     % Verasonics TX object, defining transmit parameters
        
        % Some helpful parameters
        Fs                     % The sampling frequency in Hz
        f0                     % The center frequency in Hz
        c0                     % The speed of sound in m/s
        lambda                 % The wavelength in m
        
        % For CPW and pha
        angles                 % The transmit angles of the planewaves in radians, or for pha the transmit beam angles
    end
    
    %% Constructor
    methods (Access = public)
        function h = verasonics()
            %empty constructor
        end
    end
    
    %% Set methods
    
    methods 
        function set.Trans(h,Trans)
            assert(strcmp(Trans.units,'mm'),'Please use mm as units in Verasonics.');
            h.Trans = Trans;
            h.f0 = Trans.frequency*10^6;
        end
       
        function set.Receive(h,Receive)
            assert(isempty(h.Trans)==0,'Please set the Trans variable first.');
            h.Receive = Receive;
            h.Fs = h.f0*Receive(1).samplesPerWave;
            if isfield(Receive,'aperture') == 0 % Then this is a no-mux probe and we set this to one
               h.Receive(1).aperture = 1; 
            end
        end
        
        function set.Resource(h,Resource)
            assert(isempty(h.Trans)==0,'Please set the Trans variable first.');
            h.Resource = Resource;
            h.c0 = Resource.Parameters.speedOfSound
            h.lambda = h.c0/h.f0;
        end
    end
    
    %% Read parameters
    methods (Access = public)
        function dataset = create_USTB_dataset(h,dataset)
            if ismember(class(dataset),'cpw')
                dataset = create_cpw_dataset(h,dataset);
            elseif ismember(class(dataset),'sta')
                dataset = create_sta_dataset(h,dataset);
            elseif ismember(class(dataset),'pha')
                dataset = create_pha_dataset(h,dataset);
            else
                error('Only CPW dataset is implemented for Alpinion');
            end
        end
    end
    
    methods (Access = private)
        %%%%
        %   Save to Coherent Planewave Dataset
        %
        function cpw_dataset = create_cpw_dataset(h,cpw_dataset)
            cpw_dataset.name = ['CPW Verasonics ',h.Trans.name,' Probe'];
            cpw_dataset.format = E.signal_format.RF;
            cpw_dataset.center_frequency = double(h.Trans.frequency*10^6); %center frequency in Hz
            cpw_dataset.angle = h.angles';
            cpw_dataset.geom = h.Trans.ElementPos([1:128]+h.Receive(1).aperture-1,1:3)*1e-3;
            cpw_dataset.c0 = h.c0;
            
            no_samples = h.Receive(1).endSample;
            data = zeros(no_samples, h.Resource.Parameters.numRcvChannels, length(cpw_dataset.angle), h.Resource.RcvBuffer(1).numFrames);
            delay_x0=sqrt(sum((cpw_dataset.geom-ones(128,1)*[0 0 20e-3]).^2,2))/h.c0;
            
            offset_time = calculate_delay_offset(h);
            
            % convert data
            n=1;
            cpw_dataset.time=[0:(1/h.Fs):((no_samples-1)/h.Fs)]';
            plot_delayed_signal=0;
            for n_frame = 1:h.Resource.RcvBuffer(1).numFrames
                for n_tx = 1:length(h.angles)
                    % compute time vector for this line
                    t_ini=2*h.Receive(n).startDepth*h.lambda/h.c0;
                    t_end=2*h.Receive(n).endDepth*h.lambda/h.c0;
                    no_t=(h.Receive(n).endSample-h.Receive(n).startSample+1);
                    
                    % Find t_0
                    if 1  %Calculate geometrically
                        D = abs(h.Trans.ElementPos(1,1)-h.Trans.ElementPos(end,1))*1e-3;
                        q = abs((D/2)*sin(cpw_dataset.angle(n_tx)));
                        t0_1 = q/(cpw_dataset.c0);
                    else  %Calculate using Verasonics transmit delay, this will not work for the multiplexer probe NBNB!
                        t0_1 = mean(h.TX(n_tx).Delay)*h.lambda/h.Resource.Parameters.speedOfSound;
                    end
                    
                    t_in=linspace(t_ini,t_end,no_t)-offset_time-t0_1;
                    
                    if isfield(h.Trans,'HVMux') % If Transducer has MUX, for example the L12-4v, we need to re arrange channels
                        validChannels = h.Trans.HVMux.Aperture(:,h.Receive(n_tx).aperture)';
                        validChannels = validChannels(validChannels>0);
                    else
                        validChannels = [1:128];
                    end
                    % read data
                    data(:,:,n_tx,n_frame)=interp1(t_in,double(h.RcvData{1}(h.Receive(n).startSample:h.Receive(n).endSample,validChannels,n_frame)),cpw_dataset.time,'linear',0);
                    n=n+1;

                    % to check delay calculation
                    if plot_delayed_signal
                        delay= 20e-3*cos(h.angles(n_tx))/h.c0+delay_x0;
                        
                        figure(101); hold off;
                        pcolor(1:length(cpw_dataset.geom),cpw_dataset.time,abs(data(:,:,n_tx,n_frame))); shading flat; colormap gray; colorbar; hold on;
                        plot(1:length(cpw_dataset.geom),delay,'r');
                        title(n_tx);
                        ylim([0.95*min(delay) 1.05*max(delay)]);
                        pause();
                    end
                end
            end
            
            cpw_dataset.data = data;
            
        end
        
        
        %%%%
        %    Save to Synthetic Transmit Aperture dataset
        %
        function sta_dataset = create_sta_dataset(h,sta_dataset)
            sta_dataset.name = ['STA Verasonics ',h.Trans.name,' Probe'];
            sta_dataset.format = E.signal_format.RF;
            sta_dataset.center_frequency = double(h.Trans.frequency*10^6); %center frequency in Hz
            sta_dataset.geom = h.Trans.ElementPos(1:128,1:3)*1e-3;
            sta_dataset.c0 = h.c0;
           
            no_samples=h.Receive(1).endSample;
            data=zeros(no_samples, 128, 128, h.Resource.RcvBuffer(1).numFrames);
            
            delay_x0=sqrt(sum((sta_dataset.geom-ones(128,1)*[0 0 20e-3]).^2,2))/h.c0;
            no_Apert = size(h.TX,2);
            
            offset_time = calculate_delay_offset(h);
            
            % convert data
            sta_dataset.time=[0:(1/h.Fs):((no_samples-1)/h.Fs)]';
            plot_delayed_signal=1;
            n=1;
            for n_frame = 1:h.Resource.RcvBuffer(1).numFrames
                for n_tx = 1:no_Apert
                    % compute time vector for this line
                    t_ini=2*h.Receive(n).startDepth*h.lambda/h.c0;
                    t_end=2*h.Receive(n).endDepth*h.lambda/h.c0;
                    no_t=(h.Receive(n).endSample-h.Receive(n).startSample+1);
                    t_in=linspace(t_ini,t_end,no_t)-offset_time;%scanOffsetTime;
                    
                    
                    
                    % read data
                    data(:,:,n_tx,n_frame)=interp1(t_in,double(h.RcvData{1}(h.Receive(n).startSample:h.Receive(n).endSample,:,n_frame)),sta_dataset.time,'linear',0);
                    n=n+1;
                    
                    % to check delay calculation
                    if plot_delayed_signal
                        delay = delay_x0+delay_x0(n_tx);
                        
                        figure(101); hold off;
                        pcolor(1:no_Apert,sta_dataset.time,real(data(:,:,n_tx,n_frame))); shading flat; colormap gray; colorbar; hold on;
                        plot(1:no_Apert,delay,'r');
                        title(n_tx);
                        ylim([0.9*min(delay) 1.1*max(delay)]);
                        pause();
                    end
                end
            end
            
            sta_dataset.data = data;
            
        end
        
        %% %%
        %    Save to PHased Array dataset
        %
        function pha_dataset = create_pha_dataset(h,pha_dataset)
            pha_dataset.name = ['PHA Verasonics ',h.Trans.name,' Probe'];
            pha_dataset.format = E.signal_format.RF;
            pha_dataset.center_frequency = double(h.Trans.frequency*10^6); %center frequency in Hz
            pha_dataset.geom = h.Trans.ElementPos(1:h.Trans.numelements,1:3)*1e-3;
            pha_dataset.c0 = h.c0;
            pha_dataset.angle = h.angles';
            no_samples=h.Receive(1).endSample;
            data=zeros(no_samples, h.Trans.numelements, h.Resource.RcvBuffer(1).numFrames);
            
            %delay_x0=sqrt(sum((pha_dataset.geom-ones(h.Trans.numelements,1)*[0 0 20e-3]).^2,2))/h.c0;
            %no_Apert = size(h.TX,2);
            
            offset_time = calculate_delay_offset(h);
            
            % convert data
            pha_dataset.time=[0:(1/h.Fs):((no_samples-1)/h.Fs)]';
            plot_delayed_signal=0;
            n=1;
            for n_frame = 1:h.Resource.RcvBuffer(1).numFrames
                for n_tx = 1:size(h.TX,2)
                    % compute time vector for this line
                    t_ini=2*h.Receive(n).startDepth*h.lambda/h.c0;
                    t_end=2*h.Receive(n).endDepth*h.lambda/h.c0;
                    no_t=(h.Receive(n).endSample-h.Receive(n).startSample+1);
                    
                    t0_1 = mean(h.TX(n_tx).Delay)*h.lambda/h.Resource.Parameters.speedOfSound; 
                    %t0_1 = h.TX(n_tx).Delay(32)*h.Receive(1).samplesPerWave/h.Fs;
                    
                    t_in=linspace(t_ini,t_end,no_t)-offset_time-t0_1;%scanOffsetTime;
                    
                    % read data
                    data(:,:,n_tx,n_frame)=interp1(t_in,double(h.RcvData{1}(h.Receive(n).startSample:h.Receive(n).endSample,h.Trans.Connector,n_frame)),pha_dataset.time,'linear',0);
                    n=n+1;
                    
                    % to check delay calculation
                    if plot_delayed_signal
%                         delay = delay_x0+delay_x0(n_tx);
%                         
%                         figure(101); hold off;
%                         pcolor(1:no_Apert,pha_dataset.time,real(data(:,:,n_tx,n_frame))); shading flat; colormap gray; colorbar; hold on;
%                         plot(1:no_Apert,delay,'r');
%                         title(n_tx);
%                         ylim([0.9*min(delay) 1.1*max(delay)]);
%                         pause();
                    end
                end
            end
           % data(1:20,:,:,:,:)=0;
            pha_dataset.data = data;
            
        end
        
        function offset_time = calculate_delay_offset(h)
            % offset calculation
            offset_distance=(h.TW.peak)*h.lambda;   % in [m]
            if strcmp(h.Trans.units,'mm')
                offset_distance=offset_distance+2*h.Trans.lensCorrection*1e-3;
            elseif strcmp(h.Trans.units,'wavelengths')
                offset_distance=offset_distance+2*h.Trans.lensCorrection*h.lambda;
            end
            offset_time=offset_distance/h.c0;   % in [s]
        end
        
    end
end
