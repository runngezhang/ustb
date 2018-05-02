%% Reading FI data from an UFF file recorded from a Verasonics Scanner
%
% In this example we show how to read channel data from a
% UFF (Ultrasound File Format) file recorded with a Verasonics scanner.
% You will need an internet connection to download data.
%
% _by Ole Marius Hoel Rindal <olemarius@olemarius.net>
%   and Alfonso Rodriguez-Molares <alfonso.r.molares@ntnu.no>
%
%   $Last updated: 2017/10/06$

%% Checking the file is in the path
%
% To read data from a UFF file the first we need is, you guessed it, a UFF
% file. We check if it is on the current path and download it from the USTB
% websever.

clear all; close all;

% data location
url='http://ustb.no/datasets/';      % if not found downloaded from here
%filename='L7_FI_Verasonics.uff';
filename='L7_FI_155020.uff';

% checks if the data is in your data path, and downloads it otherwise.
% The defaults data path is under USTB's folder, but you can change this
% by setting an environment variable with setenv(DATA_PATH,'the_path_you_want_to_use');
tools.download(filename, url, data_path);   

%% Reading data
%
% Let's first check if we are lucky and the file allready contains
% beamformed_data that we can display.
display=true;
content = uff.index([data_path filesep filename],'/',display);


%% Channel data
% If it doesn't have any beamformed data at least it should have some
% channel_data. So let's read that.

channel_data=uff.read_object([data_path filesep filename],'/channel_data');

%%
%
% And then do the normal routine of defining the scan,
x_axis=zeros(channel_data.N_waves,1);
for n=1:channel_data.N_waves
    x_axis(n)=channel_data.sequence(n).source.x;
end
z_axis=linspace(1e-3,62e-3,512*2).';
scan=uff.linear_scan('x_axis',x_axis,'z_axis',z_axis);

%%
%
% setting up and running the pipeline
mid=midprocess.das();
mid.dimension = dimension.both();

mid.channel_data=channel_data;
mid.scan=scan;

mid.transmit_apodization.window=uff.window.scanline;

mid.receive_apodization.window=uff.window.none;
mid.receive_apodization.f_number=1.7;

b_data=mid.go();

%% Display image
%
% And finally display the image.
b_data.plot([],'Beamformed image');


%% Beamforming with MLA's
MLA = 4;

scan_MLA=uff.linear_scan('x_axis',linspace(x_axis(1),x_axis(end),length(x_axis)*MLA)','z_axis',z_axis);
%%
mid_MLA=midprocess.das();
mid_MLA.dimension = dimension.both();

mid_MLA.channel_data=channel_data;
mid_MLA.scan=scan_MLA;

mid_MLA.transmit_apodization.window=uff.window.scanline; 
mid_MLA.transmit_apodization.MLA = MLA;
mid_MLA.transmit_apodization.MLA_overlap = MLA;

mid_MLA.receive_apodization.window=uff.window.boxcar;
mid_MLA.receive_apodization.f_number=1.7;
b_data_MLA=mid_MLA.go();

%%
b_data_MLA.plot(767,'Beamformed image MLA');
ax(1) = gca;

%%
mid_MLA_unified_fix =midprocess.das();
mid_MLA_unified_fix.dimension = dimension.both();

mid_MLA_unified_fix.channel_data=channel_data;
mid_MLA_unified_fix.scan=scan_MLA;
mid_MLA_unified_fix.use_unified_fix = 1;

mid_MLA_unified_fix.transmit_apodization.window=uff.window.scanline; 
mid_MLA_unified_fix.transmit_apodization.MLA = MLA;
mid_MLA_unified_fix.transmit_apodization.MLA_overlap = MLA/2;

mid_MLA_unified_fix.receive_apodization.window=uff.window.boxcar;
mid_MLA_unified_fix.receive_apodization.f_number=1.7;
b_data_MLA_unified_fix=mid_MLA_unified_fix.go();

%%
b_data_MLA_unified_fix.plot(777,'Beamformed image MLA unified fix');
ax(2) = gca;


%%

mid_MLA_with_plane_fix=midprocess.das();
mid_MLA_with_plane_fix.dimension = dimension.both();
mid_MLA_with_plane_fix.use_PW_fix = 1;
mid_MLA_with_plane_fix.margin_in_m = 3/1000;

mid_MLA_with_plane_fix.channel_data=channel_data;
mid_MLA_with_plane_fix.scan=scan_MLA;

mid_MLA_with_plane_fix.transmit_apodization.window=uff.window.scanline; 
mid_MLA_with_plane_fix.transmit_apodization.MLA = MLA;
mid_MLA_with_plane_fix.transmit_apodization.MLA_overlap = MLA/2;

mid_MLA_with_plane_fix.receive_apodization.window=uff.window.boxcar;
mid_MLA_with_plane_fix.receive_apodization.f_number=1.7;
b_data_MLA_with_plane_fix=mid_MLA_with_plane_fix.go();

%%
b_data_MLA_with_plane_fix.plot(778,'Beamformed image MLA with fix');

ax(3) = gca;
linkaxes(ax);

%%

tx_delay_virtual_source = reshape(mid_MLA.tx_delay_hack,scan_MLA.N_z_axis,scan_MLA.N_x_axis,channel_data.N_waves);
tx_delay_unified_fix = reshape(mid_MLA_unified_fix.tx_delay_hack,scan_MLA.N_z_axis,scan_MLA.N_x_axis,channel_data.N_waves);
tx_delay_plane_fix = reshape(mid_MLA_with_plane_fix.tx_delay_hack,scan_MLA.N_z_axis,scan_MLA.N_x_axis,channel_data.N_waves);

samples_to_skip = 1;

figure(1); clf;
b_data_MLA.plot(subplot(2,3,1),'1a : Virtual source');
ax(1) = gca;
b_data_MLA_unified_fix.plot(subplot(2,3,2),'1b : Model from [2]');
ax(2) = gca;
b_data_MLA_with_plane_fix.plot(subplot(2,3,3),'1c : Virtual source + PW model');
ax(3) = gca;
linkaxes(ax);
subplot(2,3,4); imagesc(scan_MLA.x_axis(samples_to_skip:end-samples_to_skip)*1000, scan_MLA.z_axis*1000, tx_delay_virtual_source(:,samples_to_skip:end-samples_to_skip,channel_data.N_waves/2));
title('1d: Tx delay virtual source');xlabel('x [mm]');ylabel('z [mm]');colorbar; set(gca,'fontsize',14); 
subplot(2,3,5); imagesc(scan_MLA.x_axis(samples_to_skip:end-samples_to_skip)*1000, scan_MLA.z_axis*1000, tx_delay_unified_fix(:,samples_to_skip:end-samples_to_skip,channel_data.N_waves/2)); 
title('1e: Tx delay model from [2]');xlabel('x [mm]');ylabel('z [mm]');colorbar; set(gca,'fontsize',14);
subplot(2,3,6); imagesc(scan_MLA.x_axis(samples_to_skip:end-samples_to_skip)*1000, scan_MLA.z_axis*1000, tx_delay_plane_fix(:,samples_to_skip:end-samples_to_skip,channel_data.N_waves/2));
title('1f: Tx delay virtual source + PW model');xlabel('x [mm]');ylabel('z [mm]');colorbar; set(gca,'fontsize',14);
%% Lastly and example with RTB

mid_RTB=midprocess.das();
mid_RTB.dimension = dimension.receive();
mid_RTB.use_unified_fix = 1;

mid_RTB.channel_data=channel_data;
mid_RTB.scan=scan_MLA;

mid_RTB.transmit_apodization.window=uff.window.tukey50; 
mid_RTB.transmit_apodization.f_number=3;
mid_RTB.transmit_apodization.MLA = MLA;
mid_RTB.transmit_apodization.MLA_overlap = MLA;
mid_RTB.transmit_apodization.probe = channel_data.probe;
mid_RTB.transmit_apodization.sequence = channel_data.sequence;
mid_RTB.transmit_apodization.focus = scan_MLA;
tx_apod = mid_RTB.transmit_apodization.data;

mid_RTB.receive_apodization.window=uff.window.tukey25;
mid_RTB.receive_apodization.f_number=1.7;
b_data_RTB=mid_RTB.go();
b_data_RTB.plot();


%%
tx_delay = reshape(tx_apod,scan_MLA.N_z_axis,scan_MLA.N_x_axis,channel_data.N_waves);
          
figure;
subplot(211);
imagesc(tx_delay(:,:,1));

%%
compound = postprocess.coherent_compounding()
compound.input = b_data_RTB;
b_data_RTB_img = compound.go()

%%
weighting = 1./sum(tx_apod,2);

b_data_RTB_compensated = uff.beamformed_data(b_data_RTB_img);
b_data_RTB_compensated.data = b_data_RTB_img.data .* weighting;
b_data_RTB_compensated.plot();