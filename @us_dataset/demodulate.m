function demodulate(h,plot_on,modulation_frequency,bandpass_frequency_vector,downsample_frequency, implementation)
% DEMODULATE    Demodulates RF data of signal_dataset class
%
%   Syntax:
%   demodulate(plot_on,modulation_frequency,bandpass_frequency_vector,downsample_frequency,demodulation_type)
%       plot_on                     Flag to plot signal power spectrum at all steps (Optional)
%       modulation_frequency        Modulation frequency in Hz. (Optional)
%       bandpass_frequency_vector   Vector with the bandpass frequencies in Hz [low_f_off low_f_on high_f_on high_f_off] (optional)
%       downsample_frequency        Sample frequency of the demodulated signal in Hz. (Optional) 
%       demodulation_type           Enumeration type indicating the implementation code
%
%   See also DATASET
%
% date:     01.10.2015
% authors:  Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%           Øyvind Krøvel-Velle Standal <oyvind.standal@ntnu.no>
    
    assert(h.format==E.signal_format.RF,'Signal format is not RF. Only RF data can be demodulated.');
    assert(~isempty(h.data),'Dataset is empty');
    assert(any(h.data(:)>0),'Dataset is zero');
    
    if ~exist('plot_on') plot_on=false; end
    if ~exist('implementation') implementation=E.demodulation_type.original; end
    
    % common
    fs=1/(h.time(2)-h.time(1));
                
    switch implementation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ORIGINAL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case E.demodulation_type.original
            % 2015-05-15 Alfonso Rodriguez-Molares
            
            %% power spectrum
            [fx pw] = tools.power_spectrum(h.data,h.sampling_frequency);
            assert(sum(pw)>0,'Dataset is zero');
            
            %% computing central frequency and bandwidth
            disp('Estimating power spectrum');
            fpw=filter(ones(1,26)./26,1,pw);fpw=[fpw(13:end); zeros(12,1)]; 
            [dc ic]=max(fpw.*(fx>0).'); fc=fx(ic); 
            bw_up=min(fx((fx>fc)&(fpw<dc/2).')); % -6dB upper limit
            bw_do=max(fx((fx<fc)&(fpw<dc/2).')); % -6dB down limit
            fc=(bw_up+bw_do)/2;                  % center frequency
            bw=2*(bw_up-fc);

            if ~exist('modulation_frequency') modulation_frequency=fc; end

            if(plot_on)
                assert(sum(pw)>0,'Dataset is zero');
                figure();
                subplot(1,2,1);
                plot(fx*1e-6,pw,'k'); hold on; axis manual; grid on;
                plot([modulation_frequency modulation_frequency]*1e-6,[0 1],'r--');
                plot(-[modulation_frequency modulation_frequency]*1e-6,[0 1],'r--');
                axis([-2*bw_up*1e-6 2*bw_up*1e-6 0 1]);
                xlabel('f [MHz]');
                ylabel('Relative amplitude');
                title('Before demodulation');
            end

            %% band pass filtering
            disp('Band Pass filtering');
            if ~exist('bandpass_frequency_vector') 
                transition=bw/10;
                low_freq=max([0 fc-2*bw]);      
                high_freq=min([h.sampling_frequency/2*0.99 fc+2*bw]);
                bandpass_frequency_vector=[low_freq low_freq+transition high_freq-transition high_freq];
            end
            [data fre_res w]= tools.band_pass(h.data,h.sampling_frequency,bandpass_frequency_vector);

            if(plot_on)
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                plot(fs*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
                plot(fx,pw,'r--');
                plot([modulation_frequency modulation_frequency],[0 1],'m:');
            end

            %% demodulation
            mod_sig=exp(-j*2*pi*modulation_frequency*h.time)*ones(1,size(data,2)); % demodulation sognal
            wb = waitbar(0, 'Demodulating...');
            for f=1:size(data,4)
                for n=1:size(data,3)
                    waitbar(n / size(data,3), wb);
                    data(:,:,n,f)=data(:,:,n,f).*mod_sig;
                end
            end
            close(wb);

            if(plot_on)    
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                subplot(1,2,2);
                plot(fx*1e-6,pw,'k-'); hold on; grid on; axis manual;
                axis([-2*bw_up*1e-6 2*bw_up*1e-6 0 1]);
                xlabel('f [MHz]');
                ylabel('Relative amplitude');
                title('After demodulation');
            end

            %% low pass filtering
            disp('Base Band filtering');
            baseband_frequency_vector=[0.9*modulation_frequency modulation_frequency];
            [data fre_res w]=tools.low_pass(data,h.sampling_frequency,baseband_frequency_vector);

            if(plot_on)    
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                subplot(1,2,2);
                plot(fx*1e-6,pw,'g--'); hold on; axis manual; grid on;
                plot(fs*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
                plot(-fs*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
            end

            %% resampling
            if ~exist('downsample_frequency')
                downsample_frequency=round(4*bw/1e6)*1e6; % new sampling frequency
            end

            dt=1/downsample_frequency;
            t=(h.time(1):dt:h.time(end));
            data=reshape(data,size(h.data)); % we get back singleton dimmensions
            downsample_data=zeros(length(t),size(data,2),size(data,3),size(data,4));
            wb = waitbar(0, 'Resampling...');
            for f=1:size(data,4)
                for ntx=1:size(data,3)
                    for nrx=1:size(data,2)
                        downsample_data(:,nrx,ntx,f)=interp1(h.time,data(:,nrx,ntx,f),t,'linear',0);
                    end
                    waitbar(ntx / size(data,3), wb);
                end
            end
            close(wb);

            if(plot_on)  
                [fx pw] = tools.power_spectrum(downsample_data,downsample_frequency);
                plot(fx*1e-6,pw,'m:'); 
            end

            %% write in the output 
            h.format=E.signal_format.IQ;
            h.modulation_frequency=modulation_frequency;
            h.time=t.';
            h.data=downsample_data;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% FASTFON
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case E.demodulation_type.fastfon
            % 2015-09-14 Alfonso Rodriguez-Molares
            
            siz = size(h.data);
            N = prod(siz(3:end));        % collapse dimensions beyond 3
  
            %% copy data
            data = h.data(:,:,1:N);
            
            %% show spectrum
            if(plot_on)
                [fx pw] = tools.power_spectrum(h.data,h.sampling_frequency);
                assert(sum(pw)>0,'Dataset is zero');
                figure();
                subplot(1,2,1);
                plot(fx*1e-6,pw,'k'); hold on; axis manual; grid on;
                plot([modulation_frequency modulation_frequency]*1e-6,[0 1],'r--');
                plot(-[modulation_frequency modulation_frequency]*1e-6,[0 1],'r--');
                axis([-4*modulation_frequency*1e-6 4*modulation_frequency*1e-6 0 1]);
                xlabel('f [MHz]');
                ylabel('Relative amplitude');
                title('Before demodulation');
            end
            
            %% demodulation
            if ~exist('modulation_frequency') error('Modulation frequency must be provided'); end
            mod_sig=repmat(exp(-j*2*pi*modulation_frequency*h.time),1,size(data,2),size(data,3));
            data=data.*mod_sig;
            
            if(plot_on)    
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                subplot(1,2,2);
                plot(fx*1e-6,pw,'k'); hold on; axis manual; grid on;
                plot([0 0]*1e-6,[0 1],'r--');
                plot(-2*[modulation_frequency modulation_frequency]*1e-6,[0 1],'r--');
                axis([-4*modulation_frequency*1e-6 4*modulation_frequency*1e-6 0 1]);
                xlabel('f [MHz]');
                ylabel('Relative amplitude');
                title('After demodulation');
            end

            %% low pass filtering
            %if ~exist('bandpass_frequency_vector') error('The bandpass filter frequency vector must be provided'); end
            disp('Base Band filtering');
            baseband_frequency_vector=[0.9*modulation_frequency modulation_frequency];
            [data fre_res w]=tools.low_pass(data,h.sampling_frequency,baseband_frequency_vector);

            if(plot_on)    
                [fx pw] = tools.power_spectrum(data,h.sampling_frequency);
                subplot(1,2,2);
                plot(fx*1e-6,pw,'g--'); hold on; axis manual; grid on;
                plot(fs*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
                plot(-fs*w/2/pi*1e-6,abs(fre_res)/max(abs(fre_res)),'b-')
            end
            
            %% resampling
            if ~exist('downsample_frequency') error('Downsampling frequency must be provided'); end
            Ndown=round(fs/downsample_frequency);
            
            ind_new = 1:Ndown:siz(1);       % decimating vector
            L = length(ind_new);            % length decimating  vector       
            h.time=h.time(ind_new);                             % resampled time vector
            h.data=reshape(data(ind_new,:,:),[L siz(2:end)]);   % resampled data
            h.format=E.signal_format.IQ;
            h.modulation_frequency=modulation_frequency;
                        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ØYVIND'S
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        case E.demodulation_type.oyvind
            % 2015-03-26 Øyvind K.-V. Standal
            % 2015-06-04 Alfonso Rodriguez-Molares

            doresample = false; % need to interpolate (not just downsample by integer rate)
            if ~exist('downsample_frequency'),
                Ndown = 1; % no downsampling
            else
                Ndown = round(fs/downsample_frequency);
                %if abs(Ndown - round(Ndown)) < 0.01,
                %    Ndown = round(Ndown);
                %else
                %    doresample = true;
                %    warning('Downsampling by a non-integer rate is highly inefficient.');
                %end
            end

            siz = size(h.data);
            N = prod(siz(3:end)); % collapse dimensions beyond 3
            ind_new = 1:Ndown:siz(1);
            L = length(ind_new); % length of downsampled array

            % preallocate final array
            iq = zeros([L,siz(2:end)], 'like', 1i);

            % demodulation mixing vector
            mix = exp(-2i*pi*modulation_frequency/fs*(0:siz(1)-1)');

            % otherwise not really need for low-pass filtering
            dofilter = (Ndown > 1);
            if dofilter, % make low-pass filter
                % calculate window parameters
                [M, Wn, beta, typ] = kaiserord([0.9,1]*modulation_frequency, [1,0], [1e-2,1e-2], fs);
                b = fir1(M, Wn, typ, kaiser(M+1,beta), 'noscale'); % filter design
                assert(size(h.data,1)>3*(length(b)-1),'Too many low-pass filter coefficients! Try relaxing the filter design, or increasing the samples in the dataset.');
            end

            % preallocate temp array for Hilbert
            % H = zeros(siz(1:2), 'like', 1i);

            for nn = 1:N, % loop over higher dimensions
                H = bsxfun(@times, hilbert(h.data(:,:,nn)), mix);
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

            %% write in the output 
            h.format=E.signal_format.IQ;
            h.modulation_frequency=modulation_frequency;
            dt=Ndown/fs;
            t=(h.time(1):dt:h.time(end));
            h.time=t.';
            h.data=iq;
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% CHEAP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        case E.demodulation_type.cheap
            % 2015-09-14 Alfonso Rodriguez-Molares
            
            warning('The demodulation frequency in the "cheap" demodulation is Fs/4');
            if exist('modulation_frequency')
                assert(modulation_frequency==fs/4,'The demodulation frequency in the "cheap" demodulation must be Fs/4');
            end
            
            siz = size(h.data);
            N = prod(siz(3:end));        % collapse dimensions beyond 3
  
            %% copy data
            data = h.data(:,:,1:N);
            
            %% show spectrum
            if(plot_on)
                fft_data = fftshift(fft(data(:,:), 1024));
                [fx pw] = tools.power_spectrum(h.data,h.sampling_frequency);
                assert(sum(pw)>0,'Dataset is zero');
                figure();
                subplot(1,2,1);
                plot(fx*1e-6,pw,'k'); hold on; axis manual; grid on;
                plot([modulation_frequency modulation_frequency]*1e-6,[0 1],'r--');
                plot(-[modulation_frequency modulation_frequency]*1e-6,[0 1],'r--');
                axis([-4*modulation_frequency*1e-6 4*modulation_frequency*1e-6 0 1]);
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
            
            if(plot_on)    
                [fx pw] = tools.power_spectrum(demod_data,h.sampling_frequency/2);
                subplot(1,2,2);
                plot(fx*1e-6,pw,'g--'); hold on; axis manual; grid on;
                plot([0 0]*1e-6,[0 1],'r--');
                plot(-2*[modulation_frequency modulation_frequency]*1e-6,[0 1],'r--');
                axis([-4*modulation_frequency*1e-6 4*modulation_frequency*1e-6 0 1]);
                xlabel('f [MHz]');
                ylabel('Relative amplitude');
                title('After demodulation');
            end

            h.time=h.time(1:4:L);                                     % resampled time vector
            h.data=reshape(demod_data,[L/4 siz(2:end)]);   % resampled data
            h.format=E.signal_format.IQ;
            h.modulation_frequency=fs/4;
                
    end

end








