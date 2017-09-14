classdef complex_demodulation_Denarie < preprocess
    %COMPLEX_DEMODULATION   Matlab implementation of complex demodulation
    %
    %   authors: Bastien Denaire
    %
    %   $Last updated: 2017/09/10$
    
    %% constructor
    methods (Access = public)
        function h=complex_demodulation_Denarie()
            h.name='Complex demodulation MATLAB';
            h.reference='www.ustb.no';
            h.implemented_by={'Bastien Denaire'};
            h.version='v1.0.5';
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
                warning('The modulation frequency is not specified. The estimated central frequency will be used.');
                
                % power spectrum
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
            
            % demodulation
            start_time = h.input.time(1);
            size_in    = size(h.input.data);
            
            % Find complex envelope (demodulate pre-env to zero)
            data_out      = hilbert(h.input.data(:,:)); % Find pre-envelope
            demodVect     = exp(-1i*2*pi*h.modulation_frequency*(start_time + (0:size(data_out,1)-1)*1/h.sampling_frequency)).';
            data_out      = reshape(bsxfun(@times, data_out, demodVect), size_in);
            
            % downsampling
            if isempty(h.downsample_frequency)
                warning('The downsampling frequency is not specified. Using 4*modulation_frequency');
                h.downsample_frequency=4*h.modulation_frequency;
            end
            t_new         = h.input.initial_time : 1/h.downsample_frequency : h.input.time(end);
            siz = size(data_out);
            data_out      = interp1(h.input.time, data_out(:,:), t_new, 'linear',0);
            data_out      = reshape(data_out, [size(data_out,1), siz(2:end)]);
            
            % copy results
            h.output=uff.channel_data(h.input);
            h.output.modulation_frequency=h.modulation_frequency;
            h.output.initial_time = t_new(1);
            h.output.sampling_frequency = 1/(t_new(2)-t_new(1));
            h.output.data=data_out;
            
            % pass reference
            output = h.output;
            
            % update hash
            h.save_hash();
        end
    end
    
    %% set methods
    methods
        %         function h=set.input(h,in_data)
        %             assert(isa(in_data,'uff.channel_data'), 'The input channel_data is not a CHANNEL_DATA class. Check HELP CHANNEL_DATA.');
        %             assert(in_data.modulation_frequency==0,sprintf('The input channel_data is already demodulated with %0.2 MHz',in_data.modulation_frequency/1e6));
        %             assert(~isempty(in_data.data),'The input channel_data is empty');
        %             assert(any(in_data.data(:)>0),'The input channel_data is zero');
        %
        %             timediff=in_data.time(2)-in_data.time(1);
        %             assert(timediff>0,'The time interval of the time vector in the channel_data is zero');
        %
        %             h.sampling_frequency=1/timediff;
        %             h.input=in_data;
        %         end
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
