function inter_dataset=delay_matlab(h)
% MATLAB implementation of the delay step of general beamformer

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%            Ole Marius Hoel Rindal (olemarius@olemarius.net)
%
%   $Last updated: 2017/04/24$


name='USTB Delaying step of GBeamformer MATLAB'; % beamformer name
version='v1.0.5';                      % & version

% checking aditional requirements

% modulation frequency
w0=2*pi*h.channel_data.modulation_frequency;

%% beamforming
tools.workbar();
% wave loop
for n_wave=1:numel(h.channel_data.sequence)
    
    % support multiple or single scans with the same code
    if numel(h.scan)==1
        current_scan=h.scan;
    else
        current_scan=h.scan(n_wave);
    end
    
    % precalculate receive apodization
    h.receive_apodization.probe=h.channel_data.probe;
    h.receive_apodization.scan=current_scan;
    rx_apo=h.receive_apodization.data;
    rx_propagation_distance=h.receive_apodization.propagation_distance;
    
    % precalculate transmit apodization according to 10.1109/TUFFC.2015.007183
    % compute lateral distance (assuming flat apertures, not accurate for curvilinear probes)
    h.transmit_apodization.sequence=h.channel_data.sequence(n_wave);
    h.transmit_apodization.scan=current_scan;
    tx_apo=h.transmit_apodization.data;
    
    % transmit delay
    if ~isinf(h.channel_data.sequence(n_wave).source.distance)
        % point sources
        TF=(-1).^(current_scan.z<h.channel_data.sequence(n_wave).source.z).*sqrt((h.channel_data.sequence(n_wave).source.x-current_scan.x).^2+(h.channel_data.sequence(n_wave).source.y-current_scan.y).^2+(h.channel_data.sequence(n_wave).source.z-current_scan.z).^2);
        % add distance from source to origin
        TF=TF+sign(cos(h.channel_data.sequence(n_wave).source.azimuth)).*h.channel_data.sequence(n_wave).source.distance;
    else
        % plane waves
        TF=current_scan.z*cos(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+current_scan.x*sin(h.channel_data.sequence(n_wave).source.azimuth)*cos(h.channel_data.sequence(n_wave).source.elevation)+h.scan.y*sin(h.channel_data.sequence(n_wave).source.elevation);
    end
    
    % receive loop
    for nrx=1:h.channel_data.N_elements
        tools.workbar(((n_wave-1)*h.channel_data.N_elements+nrx)/numel(h.channel_data.sequence)/h.channel_data.N_elements,sprintf('%s (%s)',name,version),'USTB');
        
        % create an intermediate beamformed data class
        inter_dataset(nrx,n_wave)=uff.beamformed_data();
        inter_dataset(nrx,n_wave).scan=current_scan;
        inter_dataset(nrx,n_wave).wave=h.channel_data.sequence(n_wave);
        inter_dataset(nrx,n_wave).data=zeros(current_scan.N_pixels,1,h.channel_data.N_frames);

        
        % receive delay
        RF=sqrt((h.channel_data.probe.x(nrx)-current_scan.x).^2+(h.channel_data.probe.y(nrx)-current_scan.y).^2+(h.channel_data.probe.z(nrx)-current_scan.z).^2);
        
        % total delay
        delay=(RF+TF)/h.channel_data.sound_speed;
        
        for n_frame=1:h.channel_data.N_frames
            
            % check whether is IQ or RF data
            if(w0>eps)
                % phase correction factor
                phase_shift=exp(1i.*w0*delay);
                % data
                data=h.channel_data.data(:,nrx,n_wave,n_frame);
            else
                % phase correction factor
                phase_shift=1;
                % data
                data=hilbert(h.channel_data.data(:,nrx,n_wave,n_frame));
            end
            
            % beamformed signal
            inter_dataset(nrx,n_wave).data(:,1,n_frame)=tx_apo.*rx_apo(:,nrx).*phase_shift.*interp1(h.channel_data.time,data,delay,'linear',0);
        end
    end
    
    % assign phase according to 2 times the receive propagation distance
    inter_dataset(n_wave).data=bsxfun(@times,inter_dataset(n_wave).data,exp(-j*w0*2*rx_propagation_distance/h.channel_data.sound_speed));
end
%close(wb);

end
