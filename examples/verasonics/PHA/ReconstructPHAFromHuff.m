recording=huff('data/20170308/P4_PHA_151704.h5','r');
recording.read();

%%
pha_dataset = recording.ultrasound_signal{1};
    
%%
% Create analytical signal or IQ if you want that :)
pha_dataset.make_analytical_signal();
%pha_dataset.demodulate(true,[],[],[],E.demodulation_algorithm.original);
%% Calculate delays
img_start_z = 0;
img_stop_z = 90/1000;
depth_samples = 512;

depth = linspace(img_start_z,img_stop_z,512);
Delays = zeros(pha_dataset.firings,size(depth,2),size(pha_dataset.geom,1)); % Buffer for the calculated delays

all_cos = cos(pha_dataset.angle);
all_sin = sin(pha_dataset.angle);
for ray = 1:pha_dataset.firings
    for i = 1:size(depth,2)
        z = all_cos(ray)*depth(i);
        x = pha_dataset.geom(:,1);
        x_focus = all_sin(ray)*depth(i);
        Delays(ray,i,:) = ((depth(i) + sqrt(z^2+(x_focus-x).^2))./pha_dataset.c0);
    end
end

%% Do beamforming
imageQube = complex(zeros(size(Delays,2),pha_dataset.firings,size(Delays,3)));% Buffer for the "image qube"
for ray = 1:size(pha_dataset.data,3)
    for iRx = 1:size(pha_dataset.data,2)
        switch (pha_dataset.format)
            case E.signal_format.RF 
                imageQube(:,ray,iRx) = interp1(pha_dataset.time,pha_dataset.data(:,iRx,ray),Delays(ray,:,iRx),'linear',0);
            case E.signal_format.AS
                imageQube(:,ray,iRx) = interp1(pha_dataset.time,pha_dataset.data(:,iRx,ray),Delays(ray,:,iRx),'linear',0);
            case E.signal_format.IQ
                phase_shift = (exp(1i.*2*pi*pha_dataset.modulation_frequency*Delays(ray,:,iRx)));
                imageQube(:,ray,iRx) =  transpose(phase_shift.*interp1(pha_dataset.time,pha_dataset.data(:,iRx,ray),Delays(ray,:,iRx),'linear',0));
        otherwise
            error('Unknown signal format!');
        end
    end
end


%% Create beamspace image
img_beamspace = mean(((imageQube)),3);
img_beamspace = db(abs(img_beamspace./max(img_beamspace(:))));

figure
imagesc(img_beamspace)

%% Do scan conversion
[img_scan_converted, Xs, Zs] = ...
            tools.scan_convert_image(img_beamspace,pha_dataset.angle,depth,1024,1024);

figure
imagesc(Xs*1000,Zs*1000,img_scan_converted);
axis image
xlabel('x [mm]');
ylabel('z [mm]');
title('Beamformed with USTB');
colormap gray
caxis([-60 0])