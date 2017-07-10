classdef verasonics < handle
    %Verasonics  Class for reading Verasonics research scanner data to USTB
    
    %   authors:    Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %               Alfonso Rodriques-Morales <alfonso.r.molares@ntnu.no>
    %
    %   $Date: 2017/03/16$
    
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
        frame_order            % The order of the frame in RcvData
        
        % For CPW
        angles                 % The transmit angles of the planewaves in radians, or for pha the transmit beam angles
        
        % For "superframes" used for Share Wave Elastography
        frames_in_superframe
        number_of_superframes
        %
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
            h.c0 = Resource.Parameters.speedOfSound;
            h.lambda = h.c0/h.f0;
        end
        function set.TW(h,TW)
            assert(isempty(h.Trans)==0,'Please set the Trans variable first.');
            h.TW=TW;
        end
        
        function frame_order = get.frame_order(h)
            % The order of the frames in the RcvBuffer is not necesarrily
            % in order. We need to check which one was the "last frame"
            if exist('h.Resource.RcvBuffer.lastFrame')
                frame_order = [h.Resource.RcvBuffer.lastFrame+1:h.Resource.RcvBuffer.numFrames 1:h.Resource.RcvBuffer.lastFrame];
            else %If the Verasonics Vantage is ran in simulation mode we don't have the lastFrame field
                frame_order = 1:h.Resource.RcvBuffer.numFrames;
            end
        end
    end
    
    % Private methods
    methods (Access = private)
        
        % Calculate the offset time from start of transmitted pulse to
        % center of pulse and compensate for the lens correction
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
        
        % Generate a USTB probe object from the Verasonics parameters
        function prb = create_probe_object(h)
            if strcmp(h.Trans.name,'L7-4') || strcmp(h.Trans.name,'P4-2v') || strcmp(h.Trans.name,'L11-4v')
                prb=uff.linear_array();
                prb.N=h.Trans.numelements;                  % number of elements
                prb.pitch=h.Trans.spacingMm/1000;           % probe pitch in azimuth [m]
                prb.element_width=h.Trans.elementWidth/1000;   % element width [m]
            else
                error('Sorry, that probe is not supported in USTB yet.');
            end
        end
        
        function trans_delays = calculate_trans_delays(h,channel_data,n_tx)
            %% Stolen from the computeTXDelays ;)
            % Hacked to work for the Verasonics definition of linear_array transmit focus with azimuth = 0
            % Delays are returned in seconds
            angle = 0;
            FocalPt(1) = channel_data.sequence(n_tx).source.x + channel_data.sequence(n_tx).source.z * sin(angle);
            FocalPt(2) = 0.0;
            FocalPt(3) = channel_data.sequence(n_tx).source.z * cos(angle);
            % Compute distance to focal point from each active element.
            X = channel_data.probe.geometry(:,1)' - FocalPt(1);
            D = sqrt(X.*X + FocalPt(3)*FocalPt(3));
            Indices = find(logical(h.TX(n_tx).Apod));
            D = max(D) - D;
            D = D - D(Indices(end));
            
            %figure(101);
            %plot(D); hold on;
            %plot(h.TX(n_tx).Delay*h.lambda)
            
            trans_delays = D/channel_data.sound_speed;
        end
        

        
        function data_out = time_shift_data(h,data_in,t_in,t_out,interpolation_factor,channel_data)
            % First do a sinc interpolation
            t_in_interp = linspace(t_in(1),t_in(end)+((interpolation_factor-1)/interpolation_factor)*(1/channel_data.sampling_frequency),length(t_in)*interpolation_factor); 
            data_tx_interpolated = interpft(double(data_in),length(t_in)*interpolation_factor);
            %%
            %                     channel = 64;
            %                     figure(99);clf;hold all;
            %                     subplot(211);hold all
            %                     plot(t_in_interp,data_tx_interpolated(:,channel),'Displayname','interpolated');
            %                     plot(t_in,data_tx(:,channel),'Displayname','original');
            %                     subplot(212);hold all
            
            %%
            % read data
            data_out=interp1(t_in_interp,data_tx_interpolated,t_out,'linear',0);
        end
    end
end
