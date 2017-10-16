classdef fast_demodulation < preprocess
    %FAST_DEMODULATION   Faster MATLAB implementation of demodulation
    %
    %   authors: Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
    %
    %   $Last updated: 2017/09/10$
    
    %% constructor
    methods (Access = public)
        function h=demodulation()
            h.name='Fast demodulation MATLAB';
            h.reference='www.ustb.no';
            h.implemented_by={'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>'};
            h.version='v1.0.0';
        end
    end
    
    properties (Access = public)
        plot_on                     % plot intermediate graphs
        modulation_frequency        % modulation frequency [Hz]
        downsample_frequency        % sampling frequency after downsampling [Hz]
    end
    
    methods
        function output=go(h)
            
            % check if we can skip calculation
            if h.check_hash()
                output= h.output;
                return;
            end
            
            % modulation frequency
            if isempty(h.modulation_frequency)||(h.modulation_frequency<eps)
                warning('The modulation frequency is not specified. The estimated central frequency will be used. For faster demodulation please provide the modulation frequency as a parameter.');

                [fx pw] = tools.power_spectrum(h.input.data,h.sampling_frequency);
                assert(sum(pw)>0,'Dataset is zero');
                
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
            
            % collapse any dimensions beyond 3rd
            siz = size(h.input.data);
            N = prod(siz(3:end));        
  
            %% copy data
            data = h.input.data(:,:,1:N);
            
            %% show spectrum
            if(h.plot_on)
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                assert(sum(pw)>0,'Dataset is zero');
                figure();
                subplot(1,2,1);
                plot(fx*1e-6,pw,'k'); hold on; axis manual; grid on;
                plot([h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
                plot(-[h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
                axis([-4*h.modulation_frequency*1e-6 4*h.modulation_frequency*1e-6 0 1]);
                xlabel('f [MHz]');
                ylabel('Relative amplitude');
                title('Before demodulation');
            end
            
            %% demodulation
            if isempty(h.modulation_frequency) error('Modulation frequency must be provided'); end
            mod_sig=repmat(exp(-j*2*pi*h.modulation_frequency*h.input.time),1,size(data,2),size(data,3));
            data=data.*mod_sig;
            
            if(h.plot_on)    
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                subplot(1,2,2);
                plot(fx*1e-6,pw,'k'); hold on; axis manual; grid on;
                plot([0 0]*1e-6,[0 1],'r--');
                plot(-2*[h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
                axis([-4*h.modulation_frequency*1e-6 4*h.modulation_frequency*1e-6 0 1]);
                xlabel('f [MHz]');
                ylabel('Relative amplitude');
                title('After demodulation');
            end

            %% low pass filtering
            disp('Base Band filtering');
            baseband_frequency_vector=[0.9*h.modulation_frequency h.modulation_frequency];
            [data fre_res w]=tools.low_pass(data,h.sampling_frequency,baseband_frequency_vector);

            if(h.plot_on)    
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                subplot(1,2,2);
                plot(fx*1e-6,pw,'g--'); hold on; axis manual; grid on;
                plot(fs*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
                plot(-fs*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
            end
            
            %% resampling
            if isempty(h.downsample_frequency) error('Downsampling frequency must be provided'); end
            Ndown=round(h.input.sampling_frequency/h.downsample_frequency);
            
            ind_new = 1:Ndown:siz(1);       % decimating vector
            L = length(ind_new);            % length decimating  vector       
            
            % write in the output
            h.output=uff.channel_data(h.input);
            h.output.modulation_frequency=h.modulation_frequency;
            h.output.initial_time=h.input.initial_time;
            h.output.sampling_frequency=h.input.sampling_frequency/Ndown;
            h.output.data=reshape(data(ind_new,:,:),[L siz(2:end)]);
            
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
    end
    
end
