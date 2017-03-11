classdef demodulator
    %demodulator   demodulator definition
    %
    %   See also PULSE, BEAM, PHANTOM, PROBE
    
    %   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
    %   $Date: 2017/02/28 $
    
    %% public properties
    properties  (SetAccess = public)
        raw_data                    % raw data class
        plot_on                     % plot intermediate graphs
        modulation_frequency        % modulation frequency [Hz]
        bandpass_frequency_vector   % bandpass trapezoidal filter appplied to signal before beamforming [Hz]
        downsample_frequency        % sampling frequency after downsampling [Hz]
    end
    
    %% private properties
    properties  (SetAccess = private)
        sampling_frequency          % sampling frequency
    end
    
    %% constructor
    methods (Access = public)
        function h=demodulator()
            %demodulator   Constructor of demodulator class
            %
            %   Syntax:
            %   h = demodulator()
            %
            %   See also BEAM, PHANTOM, PROBE, PULSE
            
            h.plot_on= true;                                               % plot intermediate graphs
        end
    end
    
    %% user methods
    methods
        
        %% go
        function out_dataset=go(h,demodulation_handle)
            if nargin<2
                out_dataset=h.original();
            else
                out_dataset=demodulation_handle();
            end
        end
        
        %% fill the gaps
        function [h,out_dataset]=fill_gaps(h)
            
            % modulation frequency
            if isempty(h.modulation_frequency)
                warning('The modulation frequency is not specified. The estimated central frequency will be used.');
                
                % power spectrum
                [fx pw] = tools.power_spectrum(h.raw_data.data,h.sampling_frequency);
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
            
            % copy dataset
            out_dataset=huff.raw_data();
            out_dataset.copy(h.raw_data);            
        end
        
    end
    
    %% set methods
    methods
        function h=set.raw_data(h,in_data)
            assert(isa(in_data,'huff.raw_data'), 'The input raw_data is not a RAW_DATA class. Check HELP RAW_DATA.');
            assert(in_data.modulation_frequency==0,sprintf('The input raw_data is already demodulated with %0.2 MHz',in_data.modulation_frequency/1e6));
            assert(~isempty(in_data.data),'The input raw_data is empty');
            assert(any(in_data.data(:)>0),'The input raw_data is zero');
            
            timediff=in_data.time(2)-in_data.time(1);
            assert(timediff>0,'The time interval of the time vector in the raw_data is zero');
            
            h.sampling_frequency=1/timediff;
            h.raw_data=in_data;
        end
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