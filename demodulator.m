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
        demodulation_algorithm      % demodulation algorithm
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
            h.demodulation_algorithm= E.demodulation_algorithm.original;   % demodulation algorithm
        end
    end
    
    %% user methods
    methods  
        
        %% go
        function out_dataset=go(h)
            
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
            out_dataset=raw_data();
            out_dataset.copy(h.raw_data);
            
            switch h.demodulation_algorithm
                case E.demodulation_algorithm.original
                    h.original(out_dataset);
                case E.demodulation_algorithm.fastfon
                    h.fastfon(out_dataset);
                case E.demodulation_algorithm.oyvind
                    h.oyvind(out_dataset);                    
                case E.demodulation_algorithm.cheap
                    h.cheap(out_dataset);
                case E.demodulation_algorithm.fieldsim
                    h.fieldsim(out_dataset);
            end
        end
        
        %% ORIGINAL: band-pass + demodulation + low pass filter 
        function original(h,out_dataset)
            % 2015-05-15 Alfonso Rodriguez-Molares
            
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
            [data fre_res w]= tools.band_pass(h.raw_data.data,h.sampling_frequency,h.bandpass_frequency_vector);

            if(h.plot_on)
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                plot(h.sampling_frequency*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
                plot(fx,pw,'r--');
                plot([h.modulation_frequency h.modulation_frequency],[0 1],'m:');
            end

            % demodulation
            mod_sig=exp(-j*2*pi*h.modulation_frequency*h.raw_data.time)*ones(1,size(data,2)); % demodulation sognal
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
            dt=1/h.downsample_frequency;
            t=(h.raw_data.time(1):dt:h.raw_data.time(end));
            data=reshape(data,size(h.raw_data.data)); % we get back singleton dimmensions
            downsample_data=zeros(length(t),size(data,2),size(data,3),size(data,4));
            wb = waitbar(0, 'Resampling...');
            for f=1:size(data,4)
                for ntx=1:size(data,3)
                    for nrx=1:size(data,2)
                        downsample_data(:,nrx,ntx,f)=interp1(h.raw_data.time,data(:,nrx,ntx,f),t,'linear',0);
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
            out_dataset.modulation_frequency=h.modulation_frequency;
            out_dataset.initial_time=t(1);
            out_dataset.sampling_frequency=1/(t(2)-t(1));
            out_dataset.data=downsample_data;
        end
        
        %% FASTFON: Demodulation + low pass filter 
        function fastfon(h,out_dataset)
           % 2015-09-14 Alfonso Rodriguez-Molares
            
            siz = size(h.data);
            N = prod(siz(3:end));        % collapse dimensions beyond 3
  
            % copy data
            data = h.data(:,:,1:N);
            
            % show spectrum
            if(h.plot_on)
                [fx pw] = tools.power_spectrum(h.data,h.sampling_frequency);
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
            
            % demodulation
            mod_sig=repmat(exp(-j*2*pi*h.modulation_frequency*h.time),1,size(data,2),size(data,3));
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
            Ndown=round(h.sampling_frequency/h.downsample_frequency);
            ind_new = 1:Ndown:siz(1);       % decimating vector
            L = length(ind_new);            % length decimating  vector       
            
            % copy results
            out_dataset.initial_time=h.raw_data.time(ind_new(1));
            out_dataset.sampling_frequency=1/(h.raw_data.time(ind_new(2))-h.raw_data.time(ind_new(1)));
            out_dataset.data=reshape(data(ind_new,:,:),[L siz(2:end)]);   % resampled data
            out_dataset.modulation_frequency=h.modulation_frequency;
        end
        
        %% ØYVIND'S Hilbert tranform demodulation
        function oyvind(h,out_dataset)
            % 2015-03-26 Øyvind K.-V. Standal
            % 2015-06-04 Alfonso Rodriguez-Molares

            doresample = false; % need to interpolate (not just downsample by integer rate)
            if empty(h.downsample_frequency),
                Ndown = 1; % no downsampling
            else
                Ndown = round(h.sampling_frequency/h.downsample_frequency);
                %if abs(Ndown - round(Ndown)) < 0.01,
                %    Ndown = round(Ndown);
                %else
                %    doresample = true;
                %    warning('Downsampling by a non-integer rate is highly inefficient.');
                %end
            end

            siz = size(h.raw_data.data);
            N = prod(siz(3:end)); % collapse dimensions beyond 3
            ind_new = 1:Ndown:siz(1);
            L = length(ind_new); % length of downsampled array

            % preallocate final array
            iq = zeros([L,siz(2:end)], 'like', 1i);

            % demodulation mixing vector
            mix = exp(-2i*pi*h.modulation_frequency/h.sampling_frequency*(0:siz(1)-1)');

            % otherwise not really need for low-pass filtering
            dofilter = (Ndown > 1);
            if dofilter, % make low-pass filter
                % calculate window parameters
                [M, Wn, beta, typ] = kaiserord([0.9,1]*h.modulation_frequency, [1,0], [1e-2,1e-2], h.sampling_frequency);
                b = fir1(M, Wn, typ, kaiser(M+1,beta), 'noscale'); % filter design
                assert(size(h.raw_data.data,1)>3*(length(b)-1),'Too many low-pass filter coefficients! Try relaxing the filter design, or increasing the samples in the dataset.');
            end

            % preallocate temp array for Hilbert
            % H = zeros(siz(1:2), 'like', 1i);

            for nn = 1:N, % loop over higher dimensions
                H = bsxfun(@times, hilbert(h.raw_data.data(:,:,nn)), mix);
                if dofilter,
                    H(:,:) = filtfilt(b, 1, H); % lowpass filter
                end
                if doresample,
                    iq(:,:,nn) = interp1(H, ind_new, 'linear', 0);
                elseif Ndown >= 2,
                    iq(:,:,nn) = H(ind_new,:); % downsample
                else
                    iq(:,:,nn) = H;
                end
            end

            % downsampling
            dt=Ndown/h.sampling_frequency;
            t=(h.time(1):dt:h.time(end));
            
            % copy results
            out_dataset.modulation_frequency=h.modulation_frequency;
            out_dataset.time=t.';
            out_dataset.raw_data.data=iq;
        end
        
        %% CHEAP 2-samples-per_wavelength demodulation
        function cheap(h,out_dataset)
            % 2015-09-14 Alfonso Rodriguez-Molares
            
            warning('The demodulation frequency in the "cheap" demodulation is h.sampling_frequency/4');
            if ~isempty(h.modulation_frequency)
                assert(h.modulation_frequency==h.sampling_frequency/4,'The demodulation frequency in the "cheap" demodulation must be h.sampling_frequency/4');
            end
            
            siz = size(h.data);
            N = prod(siz(3:end));        % collapse dimensions beyond 3
  
            %% copy data
            data = h.data(:,:,1:N);
            
            %% show spectrum
            if(h.plot_on)
                %fft_data = fftshift(fft(data(:,:), 1024));
                [fx pw] = tools.power_spectrum(h.data,h.sampling_frequency);
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
            L=siz(1)-mod(siz(1),4);
            demod_data = zeros(size(data(1:4:L,:,:)));
            demod_data(1:L/4,:,:) =     (data(1:4:L,:,:) + data(2:4:L,:,:)/j);
            %demod_data(1:2:L/2,:,:) =     data(1:4:L,:,:) + exp(-j/5)*data(2:4:L,:,:)/j;
            %demod_data(2:2:L/2,:,:) = - ( data(3:4:L,:,:) + exp(-j/5)*data(4:4:L,:,:)/j);
            
            if(h.plot_on)    
                [fx pw] = tools.power_spectrum(demod_data,h.sampling_frequency/2);
                subplot(1,2,2);
                plot(fx*1e-6,pw,'g--'); hold on; axis manual; grid on;
                plot([0 0]*1e-6,[0 1],'r--');
                plot(-2*[h.modulation_frequency h.modulation_frequency]*1e-6,[0 1],'r--');
                axis([-4*h.modulation_frequency*1e-6 4*h.modulation_frequency*1e-6 0 1]);
                xlabel('f [MHz]');
                ylabel('Relative amplitude');
                title('After demodulation');
            end

            % copy dataset
            out_dataset.time=h.time(1:4:L);                          % resampled time vector
            out_dataset.data=reshape(demod_data,[L/4 siz(2:end)]);   % resampled data
            out_dataset.modulation_frequency=h.sampling_frequency/4;
        end
        
        %% Fieldsim Hilbert tranform demodulation
        function fieldsim(h,out_dataset)
            % demodulation
            start_time = h.raw_data.time(1);
            size_in    = size(h.raw_data.data);

            % Find complex envelope (demodulate pre-env to zero)
            data_out      = hilbert(h.raw_data.data(:,:)); % Find pre-envelope
            demodVect     = exp(-1i*2*pi*h.modulation_frequency*(start_time + (0:size(data_out,1)-1)*1/h.sampling_frequency)).';
            data_out      = reshape(bsxfun(@times, data_out, demodVect), size_in);

            % downsampling
            if isempty(h.downsample_frequency) 
               warning('The downsampling frequency is not specified. Using 4*modulation_frequency'); 
               h.downsample_frequency=4*h.modulation_frequency;
            end
            t_new         = 0 : 1/h.downsample_frequency : h.time(end);
            siz = size(data_out);
            data_out      = interp1(h.time, data_out(:,:), t_new, 'linear',0);
            data_out      = reshape(data_out, [size(data_out,1), siz(2:end)]);
        
            % update dataset
            out_dataset.time=t_new.';
            out_dataset.data=data_out;
            out_dataset.modulation_frequency=h.modulation_frequency;        
        end
    end
    
    %% set methods
    methods  
        function h=set.raw_data(h,in_data)
            assert(strcmp(class(in_data),'raw_data'), 'The input raw_data is not a RAW_DATA class. Check HELP RAW_DATA.');
            assert(in_data.modulation_frequency==0,sprintf('The input raw_data is already demodulated with %0.2 MHz',in_data.modulation_frequency/1e6));
            assert(~isempty(in_data.data),'The input raw_data is empty');
            assert(any(in_data.data(:)>0),'The input raw_data is zero');

            timediff=in_data.time(2)-in_data.time(1);
            assert(timediff>0,'The time interval of the time vector in the raw_data is zero');
            
            h.sampling_frequency=1/timediff;
            h.raw_data=in_data;
        end
        function h=set.plot_on(h,in_plot_on)
            assert(strcmp(class(in_plot_on),'logical'), 'The input plot_on is not a LOGICAL class (true/false). Check HELP LOGICAL.');
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
        function h=set.demodulation_algorithm(h,in_demodulation_algorithm)
            assert(strcmp(class(in_demodulation_algorithm),'E.demodulation_algorithm'), 'The input demodulation_algorithm is not a E.demodulation_algorithm class. Check HELP E.demodulation_algorithm.');
            h.demodulation_algorithm=in_demodulation_algorithm;
        end        
    end
end