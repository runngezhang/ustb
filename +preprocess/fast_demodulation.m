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
        baseband_frequency_vector = [0.5, 1];       % start and end of the transition band, defined
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
                fprintf(1, 'Estimating power spectrum\n');
                [dc, ic] = max(pw); 

                bw_lo = interp1(pw(1:ic), fx(1:ic), dc/2);          % -6dB lower limit
                bw_up = interp1(pw(ic:end/2), fx(ic:end/2), dc/2);  % -6dB upper limit
                fc = (bw_lo+bw_up)/2;                               % center frequency
                
                % Set modulation ferquency
                h.modulation_frequency = -fc;
            end
            
            % estimate downsampling frequency if needed
            if isempty(h.downsample_frequency)
                warning(['The downsampling frequency is not specified. ', ...
                    'Using 2 * modulation_frequency.'])
                
                h.downsample_frequency = 2 * h.modulation_frequency;
            end
            
            if mod(h.input.sampling_frequency, h.downsample_frequency) ~= 0
                warning(['The input sample frequency is not a multiple of the specified downsample frequency. ', ...
                    'Rounding up to the next higher downsample frequency.'])
            end
            
            Ndown = floor(h.input.sampling_frequency / h.downsample_frequency);
            h.downsample_frequency = h.input.sampling_frequency / Ndown;
            
            % Plot RF channel data power spectrum
            if(h.plot_on)      
                if ~exist('pw', 'var')
                    [fx, pw] = tools.power_spectrum(h.input.data, h.sampling_frequency);
                end
                
                pv = max(pw);       % find peak value
                
                figure('Color', 'w')
                subplot(1,2,1)
                hold on
                obj = plot(fx*1e-6, 10*log10(pw), 'k', 'LineWidth', 1, 'DisplayName', 'RF channel data');
                plot([h.modulation_frequency, h.modulation_frequency]*1e-6, [-120, 0]+10*log10(pv), 'r--', ...
                    'LineWidth', 1)
                plot(-[h.modulation_frequency, h.modulation_frequency]*1e-6, [-120, 0]+10*log10(pv), 'r--', ...
                    'LineWidth', 1)
                hold off
                xlim([-2*h.modulation_frequency, 2*h.modulation_frequency]*1e-6)
                ylim([-120, 0])
                grid on
                box on
                xlabel('f [MHz]')
                ylabel('Power spectrum [dB]')
                legend(obj, 'location', 'southeast')
            end
            
            % Demodulation
            data = h.input.data .* exp(-1j*2*pi*h.modulation_frequency*h.input.time);
            
            if(h.plot_on)    
                [fx, pw] = tools.power_spectrum(data, h.sampling_frequency);
                
                pv = max(pw);       % find peak value

                subplot(1,2,2)
                hold on
                obj = plot(fx*1e-6, 10*log10(pw), 'k', 'LineWidth', 1, 'DisplayName', 'Down-mixed channel data');
                plot([0, 0]*1e-6, [-120, 0]+10*log10(pv), 'r--', 'LineWidth', 1)
                hold off
                xlim([-h.downsample_frequency, h.downsample_frequency]*1e-6)
                ylim([-120, 0])
                grid on
                box on
                xlabel('f [MHz]');
                ylabel('Power spectrum [dB]');
            end

            % Perform base-band filtering
            fprintf(1, 'Base Band filtering\n');
            [data, H, W, delay] = tools.low_pass(data, h.sampling_frequency, ...
                h.baseband_frequency_vector*h.modulation_frequency);
            
            if(h.plot_on)    
                [fx, pw] = tools.power_spectrum(data, h.sampling_frequency);
                
                pv = max(pw);       % find peak value

                subplot(1,2,2)
                hold on
                obj(2) = plot(fx*1e-6, 10*log10(pw), 'c--', 'LineWidth', 1, 'DisplayName', 'Base-band filtered channel data');
                obj(3) = plot(h.input.sampling_frequency*W/2/pi*1e-6, 10*log10(abs(H).^2 / max(abs(H).^2)) + 10*log10(pv), 'b-', 'DisplayName', 'Base-band filter frequency response');
                plot(-h.input.sampling_frequency*W/2/pi*1e-6, 10*log10(abs(H).^2 / max(abs(H).^2)) + 10*log10(pv), 'b-')
                hold off
                legend(obj, 'location', 'southeast')
            end
            
            % Downsampling            
            id = delay + 1:Ndown:size(data, 1);       % decimating vector
            
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
        function set.plot_on(h, val)
            validateattributes(val, {'logical'}, {'scalar'})
            h.plot_on = val;
        end
        function set.modulation_frequency(h, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real'})
            h.modulation_frequency = val;
        end
        function set.downsample_frequency(h, val)
            validateattributes(val, {'numeric'}, {'scalar', 'real'})
            h.downsample_frequency = val;
        end
        function set.baseband_frequency_vector(h, val)
            validateattributes(val, {'numeric'}, {'nonnegative', 'size', [1, 2], 'real'})
            h.baseband_frequency_vector = val;
        end
    end
end
