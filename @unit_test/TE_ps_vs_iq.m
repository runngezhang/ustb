function ok = TE_ps_vs_iq(h)
%PS_VS_IQ Point Spread function Diverging Waves IQ test
%   Downloads data from 'http://hirse.medisin.ntnu.no/ustb/data/ps/'
%   beamforms it and compares it with previously beamformed data (USTB v1.9)

    import huff.*;
    
    % data location
    url='http://hirse.medisin.ntnu.no/ustb/data/ps/';   % if not found data will be downloaded from here
    local_path='data/ps/';                              % location of example data in this computer                      
    raw_data_filename='ps_vs_iq.mat';
    beamformed_data_filename='beamformed_ps_vs_iq.mat';
    
    % check if the file is available in the local path & downloads otherwise
    h.download(raw_data_filename, url, local_path);
    h.download(beamformed_data_filename, url, local_path);
    
    % load data
    load([local_path raw_data_filename]);    
    load([local_path beamformed_data_filename]);    
    
    % PROBE
    prb=probe(s.geom);
    
    % SEQUENCE 
    for n=1:length(s.source)
        seq(n)=wave();
        seq(n).probe=prb;
        seq(n).sound_speed=s.c0;
        seq(n).source.xyz=s.source(n,:);
    end
    
    % converting time vector standards
    data=zeros(size(s.data));
    for n=1:length(s.source)
        delay=-seq(n).source.distance/s.c0;
        pcf=exp(-1i.*2*pi*s.modulation_frequency*delay);
        data(:,:,n)=pcf.*interp1(s.time+delay,s.data(:,:,n),s.time,'pchip',0);
    end
        
    % RAW DATA
    r_data=channel_data();
    r_data.probe=prb;
    r_data.sequence=seq;
    r_data.initial_time=s.time(1);
    r_data.sampling_frequency=1/(s.time(2)-s.time(1));
    r_data.modulation_frequency=s.modulation_frequency;
    r_data.sound_speed=s.c0;
    r_data.data=data;
    
    % APODIZATION
    apo=apodization();
    apo.window=window.boxcar;
    apo.f_number=r.f_number;
    apo.apex.distance=Inf;
    
    % BEAMFORMER
    bmf=beamformer();
    bmf.channel_data=r_data;
    bmf.receive_apodization=apo;
    bmf.transmit_apodization=apo;
    bmf.scan=linear_scan(r.x_axis,r.z_axis);
        
    % beamforming
    b_data=bmf.go(@postprocess.coherent_compound);

    % test result
    ok=(norm(b_data.data-r.data(:))/norm(r.data(:)))<h.external_tolerance;

%     figure;
%     plot(real(b_data.data)); hold on; grid on;
%     plot(real(r.data(:)),'r--'); 
%    
%     % show
%     b_data.plot([],'Result');
%     
%     % ref
%     ref_data=beamformed_data();
%     ref_data.data=r.data;
%     ref_data.scan=linear_scan(r.x_axis,r.z_axis);
%     ref_data.plot([],'Reference');
    
end

