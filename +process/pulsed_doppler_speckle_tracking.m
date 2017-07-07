classdef pulsed_doppler_speckle_tracking < process
    %COHERENT_COMPOUNDING   Matlab implementation of Coherent compounding
    %
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %            Ole Marius Hoel Rindal <olemarius@olemarius.net>
    %
    %   $Last updated: 2017/05/11$
    
    %% constructor
    methods (Access = public)
        function h=pulsed_doppler_speckle_tracking()
            h.name='Pulsed Doppler Speckle Tracking';
            h.reference='www.USTB.no'; %Need to add some references here!
            h.implemented_by={'Ole Marius Hoel Rindal <olemarius@olemarius.net>'};
            h.version='v1.0.0';
        end
    end
    
    %% Additional properties
    properties
        z_gate = 4
        x_gate = 2
        packet_size = 6
    end
    
    methods
        function out_data=go(h)
            [N_pixels Nrx Ntx N_frames]=size(h.beamformed_data.data);
            
            assert(N_frames>h.packet_size,'The number of frames needs to be higher than the packet size');
            assert(Nrx==1,'The pulsed doppler speckle traking can only be used between frames');
            assert(Ntx==1,'The pulsed doppler speckle traking can only be used between frames');
            assert(mod(h.z_gate,2)==0,'Please use an even number for the z_gate');
            assert(mod(h.x_gate,2)==0,'Please use an even number for the x_gate');
            
            % declare output structure
            out_data=uff.beamformed_data(h.beamformed_data); % ToDo: instead we should copy everything but the data
            
            % calculate sampling frequency in image
            h.beamformed_data.calculate_sampling_frequency(h.channel_data.sound_speed);
            
            % get images in matrix format
            images = h.beamformed_data.get_image('none-complex');
            
            % create a buffer for the output
            displacement_data = zeros(size(h.beamformed_data.data,1),size(h.beamformed_data.data,2),...
                        size(h.beamformed_data.data,3),size(h.beamformed_data.data,4)-h.packet_size+1);
                    
            for i = h.packet_size:h.beamformed_data.N_frames
                % buffer
                temp_disp = zeros(size(images(:,:,1,1,1)));
                
                %Calculate displacement
                [temp_disp(h.z_gate/2:end-h.z_gate/2-1,h.x_gate/2:end-h.x_gate/2)] = pulsed_doppler_displacement_estimation(h,images(:,:,i-h.packet_size+1:i));
                displacement_data(:,1,1,i-h.packet_size+1) = temp_disp(:);
            end
            
            out_data.data = displacement_data;
        end
    end
    
    methods (Access = private)
        function [d_1,d_2,f_hat,fc_hat,C] = pulsed_doppler_displacement_estimation(h,X)
            
            c   = h.channel_data.sound_speed;           %Speed of sound
            fs  = h.beamformed_data.sampling_frequency; %Sampling frequency
            fc = h.channel_data.pulse.center_frequency; %Central frequency
            U = h.z_gate;
            V = h.x_gate;
            O = h.packet_size;
            
            X_conj = conj(X); %Complex conjugate of X
            
            %Autocorrelation lag
            R_0_1 = sum(X(1:end-1,:,1:end-1).*X_conj(1:end-1,:,2:end),3);
            R_0_1 = conv2(R_0_1, ones(U,1), 'valid');
            R_0_1 = conv2(R_0_1, ones(1,V), 'valid');
            
            %Autocorrelation lag
            R_1_0 = sum(X(1:end-1,:,:).*X_conj(2:end,:,:),3);
            R_1_0 = conv2(R_1_0, ones(U,1), 'valid');
            R_1_0 = conv2(R_1_0, ones(1,V), 'valid');
            
            %Estimated doppler frequency
            f_hat = (angle((R_0_1))/(2*pi));
            
            %Estimated centeral frequency
            fc_hat = abs((angle(R_1_0))/(2*pi*1/fs));
            
            %Correlation coefficient estimation qualityindicator.
            C = sum(X(1:end-1,:,:).*X_conj(1:end-1,:,:),3);
            C = conv2(C, ones(U,1), 'valid');
            C = conv2(C, ones(1,V), 'valid');
            C = (O/(O-1))*abs(R_0_1)./C;
            
            %Autocorrelation method
            d_1 = c*f_hat/(2*fc);
            %Modified autocorrelation method
            d_2 = c*f_hat./(2*fc_hat);
        end
        
        
    end
end