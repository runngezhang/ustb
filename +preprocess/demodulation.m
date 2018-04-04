classdef demodulation < preprocess
    %DEMODULATION   Matlab implementation of demodulation
    %
    %   authors: Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
    %
    %   $Last updated: 2017/09/10$
    
    %% constructor
    methods (Access = public)
        function h=demodulation()
            h.name='Demodulation MATLAB';
            h.reference='www.ustb.no';
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};
            h.version='v1.0.5';
        end
    end
    
    properties (Access = public)
        plot_on                     % plot intermediate graphs
        modulation_frequency        % modulation frequency [Hz]
        bandpass_frequency_vector   % bandpass trapezoidal filter appplied to signal before beamforming [Hz]
        downsample_frequency        % sampling frequency after downsampling [Hz]
    end
    
    methods
        function output=go(h)
            
            % check if we can skip calculation
            if h.check_hash()
                output= h.output;
                return;
            end
            
            % power spectrum
            [fx pw] = tools.power_spectrum(h.input.data,h.sampling_frequency);
            assert(sum(pw)>0,'Dataset is zero');
            
            % modulation frequency
            if isempty(h.modulation_frequency)||(h.modulation_frequency<eps)
                warning('The modulation frequency is not specified. The estimated central frequency will be used.');
                
                % computing central frequency and bandwidth
                disp('Estimating power spectrum');
                fpw=filter(ones(1,26)./26,1,pw);fpw=[fpw(13:end); zeros(12,1)];
                [dc ic]=max(fpw.*(fx>0).'); fc=fx(ic);
                bw_up=min(fx((fx>fc)&(fpw<dc/2).')); % -6dB upper limit
                bw_do=max(fx((fx<fc)&(fpw<dc/2).')); % -6dB down limit
                fc=(bw_up+bw_do)/2;                  % center frequency
                bw=2*(bw_up-fc);
                
                % set modulation ferquency
                h.modulation_frequency=fc;
            end
            
            % downsampling frequency
            if isempty(h.downsample_frequency)
                warning('The downsampling frequency is not specified. Using 4*modulation_frequency.');
                h.downsample_frequency=4*h.modulation_frequency;
            end
            
            % computing central frequency and bandwidth
            disp('Estimating power spectrum');
            fpw=filter(ones(1,26)./26,1,pw);fpw=[fpw(13:end); zeros(12,1)];
            [dc ic]=max(fpw.*(fx>0).'); fc=fx(ic);
            bw_up=min(fx((fx>fc)&(fpw<dc/2).')); % -6dB upper limit
            bw_do=max(fx((fx<fc)&(fpw<dc/2).')); % -6dB down limit
            fc=(bw_up+bw_do)/2;                  % center frequency
            bw=2*(bw_up-fc);
            
            if(h.plot_on)
                assert(sum(pw)>0,'Dataset is zero');
                figure();
                subplot(1,2,1);
                plot(fx*1e-6,pw,'k'); hold on; axis manual; grid on;
                plot([h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
                plot(-[h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
                axis([-2*bw_up*1e-6 2*bw_up*1e-6 0 1]);
                xlabel('f [MHz]');
                ylabel('Relative amplitude');
                title('Before demodulation');
            end
            
            % band pass filtering
            disp('Band Pass filtering');
            if isempty(h.bandpass_frequency_vector)
                transition=bw/10;
                low_freq=max([0 fc-2*bw]);
                high_freq=min([h.sampling_frequency/2*0.99 fc+2*bw]);
                h.bandpass_frequency_vector=[low_freq low_freq+transition high_freq-transition high_freq];
            end
            [data fre_res w]= tools.band_pass(h.input.data,h.sampling_frequency,h.bandpass_frequency_vector);
            
            if(h.plot_on)
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                plot(h.sampling_frequency*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
                plot(fx,pw,'r--');
                plot([h.modulation_frequency h.modulation_frequency],[0 1],'m:');
            end
            
            % demodulation
            mod_sig=exp(-j*2*pi*h.modulation_frequency*h.input.time)*ones(1,size(data,2)); % demodulation sognal
            wb = waitbar(0, 'Demodulating...');
            for f=1:size(data,4)
                for n=1:size(data,3)
                    waitbar(n / size(data,3), wb);
                    data(:,:,n,f)=data(:,:,n,f).*mod_sig;
                end
            end
            close(wb);
            
            if(h.plot_on)
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                subplot(1,2,2);
                plot(fx*1e-6,pw,'k-'); hold on; grid on; axis manual;
                axis([-2*bw_up*1e-6 2*bw_up*1e-6 0 1]);
                xlabel('f [MHz]');
                ylabel('Relative amplitude');
                title('After demodulation');
            end
            
            % low pass filtering
            disp('Base Band filtering');
            baseband_frequency_vector=[0.9*h.modulation_frequency h.modulation_frequency];
            [data fre_res w]=tools.low_pass(data,h.sampling_frequency,baseband_frequency_vector);
            
            if(h.plot_on)
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                subplot(1,2,2);
                plot(fx*1e-6,pw,'g--'); hold on; axis manual; grid on;
                plot(h.sampling_frequency*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
                plot(-h.sampling_frequency*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
            end
            
            % resampling
            dt=1./h.downsample_frequency;
            temp_time = h.input.time;
            t =  temp_time(1):dt:temp_time(end);
            data=reshape(data,size(h.input.data)); % we get back singleton dimmensions
            downsample_data=zeros(length(t),size(data,2),size(data,3),size(data,4));
            wb = waitbar(0, 'Resampling...');
            for f=1:size(data,4)
                for ntx=1:size(data,3)
                    for nrx=1:size(data,2)
                        downsample_data(:,nrx,ntx,f)=interp1(h.input.time,data(:,nrx,ntx,f),t,'linear',0);
                    end
                    waitbar(ntx / size(data,3), wb);
                end
            end
            close(wb);
            
            if(h.plot_on)
                [fx pw] = tools.power_spectrum(downsample_data,h.downsample_frequency);
                plot(fx*1e-6,pw,'m:');
            end
            
            % write in the output
            h.output=uff.channel_data(h.input);
            h.output.modulation_frequency=h.modulation_frequency;
            h.output.initial_time=t(1);
            h.output.sampling_frequency=1/(t(2)-t(1));
            h.output.data=downsample_data;
            
            
            % pass reference
            output = h.output;
            
            % update hash
            h.save_hash();
        end
    end
    
    %% set methods
    methods
        function h=set.plot_on(h,in_plot_on)
            assert(isa(in_plot_on,'logical'), 'The input plot_on is not a LOGICAL class (true/false). Check HELP LOGICAL.');
            h.plot_on=in_plot_on;
        end
        function h=set.modulation_frequency(h,in_modulation_frequency)
            assert(numel(in_modulation_frequency)==1, 'The modulation_frequency must be a escalar');
            h.modulation_frequency=in_modulation_frequency;
        end
        function h=set.downsample_frequency(h,in_downsample_frequency)
            assert(numel(in_downsample_frequency)==1, 'The downsample_frequency must be a escalar');
            h.downsample_frequency=in_downsample_frequency;
        end
        function h=set.bandpass_frequency_vector(h,in_bandpass_frequency_vector)
            assert(numel(in_bandpass_frequency_vector)==4, 'The bandpass_frequency_vector must have 4 elements [F_low_out F_low_in F_up_in F_up_out]');
            if size(in_bandpass_frequency_vector,1)==1
                h.bandpass_frequency_vector=in_bandpass_frequency_vector;
            else
                h.bandpass_frequency_vector=in_bandpass_frequency_vector.';
            end
        end
    end
    
end
