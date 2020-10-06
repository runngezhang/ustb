classdef fast_demodulation < preprocess
    %FAST_DEMODULATION   Faster MATLAB implementation of demodulation
    %
    %   authors: Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
    %            Stefano Fiorentini <stefano.fiorentini@ntnu.no>
    %
    %   $Last updated: 2020/10/06$
    
    %% constructor
    methods (Access = public)
        function h = fast_demodulation()
            h.name = 'Fast demodulation MATLAB';
            h.reference = 'www.ustb.no';
            h.implemented_by = {'Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>', ...
                'Stefano Fiorentini <stefano.fiorentini@ntnu.no>'};
            h.version = 'v1.0.1';
        end
    end
    
    properties (Access = public)
        plot_on                                     % plot intermediate graphs
        modulation_frequency                        % modulation frequency [Hz]
        downsample_frequency                        % sampling frequency after downsampling [Hz]
        baseband_frequency_vector = [0.5, 1.5];     % start and end of the transition band, defined
                                                    % as a multiple of the modulation frequency
    end
    
    methods
        function output=go(h)
            
            % check if we can skip calculation
            if h.check_hash()
                output= h.output;
                return;
            end
            
            % estimate modulation frequency if needed
            if isempty(h.modulation_frequency)
                warning(['The modulation frequency is not specified. ', ...
                    'The estimated central frequency will be used. '])

                [fx, pw] = tools.power_spectrum(h.input.data, h.sampling_frequency);
                
                % computing central frequency and bandwidth
                fprintf(1, 'Estimating power spectrum');
                [dc, ic] = max(pw); 

                bw_lo = interp1(pw(1:ic), fx(1:ic), dc/2);      % -6dB lower limit
                bw_up = interp1(pw(ic:end), fx(ic:end), dc/2);  % -6dB upper limit
                fc = (bw_lo+bw_up)/2;                           % center frequency
                
                % set modulation ferquency
                h.modulation_frequency = fc;
            end
            
            % estimate downsampling frequency if needed
            if isempty(h.downsample_frequency)
                warning(['The downsampling frequency is not specified. ', ...
                    'Using 2*modulation_frequency.'])
                
                h.downsample_frequency = 2*h.modulation_frequency;
            else
                if mod(h.input.sampling_frequency, h.downsample_frequency) ~= 0
                    warning(['The input sample frequency is not a multiple of the specified downsample frequency. ', ...
                        'Rounding up to the next higher downsample frequency.'])
                end
                
                Ndown = ceil(h.input.sampling_frequency / h.downsample_frequency);
                h.downsample_frequency = h.input.sampling_frequency / Ndown;
            end
              
            % Plot RF channel data power spectrum
            if(h.plot_on)      
                if ~exist('pw', 'var')
                    [fx, pw] = tools.power_spectrum(h.input.data, h.sampling_frequency);
                end
                
                pv = max(pw);       % find peak value
                
                figure('Color', 'w')
                subplot(1,2,1)
                hold on
                plot(fx*1e-6, 10*log10(pw), 'k', 'LineWidth', 1)
                plot([h.modulation_frequency, h.modulation_frequency]*1e-6, [-60, 0]+10*log10(pv), 'r--', ...
                    'LineWidth', 1)
                plot(-[h.modulation_frequency, h.modulation_frequency]*1e-6, [-60, 0]+10*log10(pv), 'r--', ...
                    'LineWidth', 1)
                hold off
                xlim([-2*h.modulation_frequency, 2*h.modulation_frequency])
                ylim([-60, 0])
                axis tight
                grid on
                box on
                xlabel('f [MHz]');
                ylabel('Power spectrum [dB]');
                title('Before demodulation');
            end
            
            % Demodulation
            data = h.input.data .* exp(-1j*2*pi*h.modulation_frequency*h.input.time);
            
            if(h.plot_on)    
                [fx, pw] = tools.power_spectrum(data, h.sampling_frequency);
                
                pv = max(pw);       % find peak value

                subplot(1,2,2)
                hold on
                plot(fx*1e-6, pw, 'k', 'LineWidth', 1)
                plot([0, 0]*1e-6, [-60, 0]+10*log10(pv), 'r--', 'LineWidth', 1)
                hold off
                xlim([-h.downsample_frequency/2, h.downsample_frequency/2])
                ylim([-60, 0])
                axis tight
                grid on
                box on
                xlabel('f [MHz]');
                ylabel('Power spectrum [dB]');
                title('Before demodulation');
            end

            % Perform base-band filtering
            fprintf(1, 'Base Band filtering');
            [data, fre_res, w, delay] = tools.low_pass(data, h.sampling_frequency, ...
                h.baseband_frequency_vector*h.modulation_frequency);

            if(h.plot_on)    
                [fx, pw] = tools.power_spectrum(data, h.sampling_frequency);
                
                pv = max(pw);       % find peak value

                subplot(1,2,2)
                hold on
                plot(fx*1e-6, 10*log10(pw), 'c--', 'LineWidth', 1)
                plot(fs*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
                plot(-fs*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
                hold on
            end
            
            % Downsampling            
            id = delay+1:Ndown:size(data, 1);       % decimating vector
            
            % Create output channel datta object
            h.output = uff.channel_data();
            h.output.modulation_frequency = h.modulation_frequency;
            h.output.initial_time = h.input.initial_time;
            h.output.sampling_frequency = h.downsample_frequency;
            h.output.data = data(id, :, :, :);
            
            % pass reference
            output = h.output;
            
            % update hash
            h.save_hash();
        end
    end
    
    %% set methods
    methods
        function set.plot_on(h, in_plot_on)
            assert(isa(in_plot_on,'logical'), 'The input plot_on is not a LOGICAL class (true/false). Check HELP LOGICAL.');
            h.plot_on = in_plot_on;
        end
        function set.modulation_frequency(h, in_modulation_frequency)
            assert(numel(in_modulation_frequency)==1, 'The modulation_frequency must be a escalar');
            h.modulation_frequency=in_modulation_frequency;
        end
        function set.downsample_frequency(h, in_downsample_frequency)
            assert(numel(in_downsample_frequency)==1, 'The downsample_frequency must be a escalar');
            h.downsample_frequency=in_downsample_frequency;
        end
    end
    
end
