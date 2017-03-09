function ok = ps_sta_iq(h)
%PS_STA_IQ Point Spread function STA IQ test
%   Downloads data from 'http://hirse.medisin.ntnu.no/ustb/data/ps/'
%   beamforms it and compares it with previously beamformed data

    % data location
    url='http://hirse.medisin.ntnu.no/ustb/data/ps/';   % if not found data will be downloaded from here
    local_path='data/ps/';                              % location of example data in this computer                      
    raw_data_filename='ps_sta_rf.mat';
    beamformed_data_filename='beamformed_ps_sta_rf.mat';
    
    % check if the file is available in the local path & downloads otherwise
    h.download(raw_data_filename, url, local_path);
    h.download(beamformed_data_filename, url, local_path);
    
    % load data
    load([local_path raw_data_filename]);    
    load([local_path beamformed_data_filename]);    
    
    % PROBE
    prb=probe(s.geom);
    
    % SEQUENCE 
    for n=1:prb.N_elements 
        seq(n)=wave();
        seq(n).probe=prb;
        seq(n).sound_speed=s.c0;
        seq(n).source.xyz=[prb.x(n) prb.y(n) prb.z(n)];
    end
    
    % converting time vector standards
    data=zeros(size(s.data));
    for n=1:prb.N_elements 
        data(:,:,n)=interp1(s.time+seq(n).source.distance/s.c0,s.data(:,:,n),s.time,'linear',0);
    end
    
    % RAW DATA
    r_data=raw_data();
    r_data.probe=prb;
    r_data.sequence=seq;
    r_data.initial_time=s.time(1);
    r_data.sampling_frequency=1/(s.time(2)-s.time(1));
    r_data.sound_speed=s.c0;
    r_data.data=data;
    
    
%     fh=figure;
%     for n=1:prb.N_elements 
%         r_data.plot(n); 
%         dst=norm(seq(n).source.xyz-[0 0 40e-3])+sqrt(sum((prb.geometry(:,1:3)-ones(prb.N_elements,1)*[0 0 40e-3]).^2,2));
%         delay=dst/s.c0+seq(n).source.distance/s.c0;
%         %subplot(1,2,1); 
%         hold on;
%         plot(delay*1e6,'r--','linewidth',2)
%         %subplot(1,2,2); hold on;
%         %plot(delay*1e6,'r--','linewidth',2)
%         pause;
%     end
    
    % APODIZATION
    apo=apodization();
    apo.window=window.boxcar;
    apo.f_number=r.f_number;
    apo.apex.distance=Inf;
    
    % BEAMFORMER
    bmf=beamformer();
    bmf.raw_data=r_data;
    bmf.receive_apodization=apo;
    bmf.transmit_apodization=apo;
    bmf.scan=linear_scan(r.x_axis,r.z_axis);
        
    % beamforming
    b_data=bmf.go(@postprocess.coherent_compound);
    
    % show
    b_data.plot([],'Result');
    
    % ref
    ref_data=beamformed_data();
    ref_data.data=r.data;
    ref_data.scan=linear_scan(r.x_axis,r.z_axis);
    ref_data.plot([],'Reference');
    
end

