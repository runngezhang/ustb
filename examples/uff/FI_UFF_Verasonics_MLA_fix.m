%% Focused Linear scan with L7-4 probe in Verasonics demonstrating MLA
%
%   This example demonstrates the MLA implementation and demonstrates
%   different fixes to the artifact occuring near the focus.
%
%   The solution with the transmit delay model introduced in  Nguyen, N. Q., & Prager, R. W. (2016). 
%   High-Resolution Ultrasound Imaging With Unified Pixel-Based Beamforming. IEEE Trans. Med. Imaging, 35(1), 98?108.
%   But also a second simpler model assuming PW around focus.
%
%
% _by Ole Marius Hoel Rindal <olemarius@olemarius.net>
%
%   $Last updated: 2018/20/04$

%% Read channeldata

clear all; close all;

% data location
url='http://ustb.no/datasets/';      % if not found downloaded from here
filename='L7_FI_Verasonics.uff';
%filename='L7_FI_155020.uff'; % The same dataset, but the foci hits the point scatteres, enhancing the artifact
tools.download(filename, url, data_path);   
channel_data=uff.read_object([data_path filesep filename],'/channel_data');

%% Define Scan
x_axis=zeros(channel_data.N_waves,1);
for n=1:channel_data.N_waves
    x_axis(n)=channel_data.sequence(n).source.x;
end
z_axis=linspace(1e-3,62e-3,512*2).';
scan=uff.linear_scan('x_axis',x_axis,'z_axis',z_axis);

%% Conventional Scanline Beamforming
mid=midprocess.das();
mid.dimension = dimension.both();

mid.channel_data=channel_data;
mid.scan=scan;

mid.transmit_apodization.window=uff.window.scanline;

mid.receive_apodization.window=uff.window.none;
mid.receive_apodization.f_number=1.7;

b_data=mid.go();

%% Display image
b_data.plot([],'Conventional one scanline per transmit');


%% Beamforming with MLA's
% create MLA scan
MLA = 4;
scan_MLA=uff.linear_scan('x_axis',linspace(x_axis(1),x_axis(end),length(x_axis)*MLA)','z_axis',z_axis);

% beamform without any fix using conventional virtual source model
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

b_data_MLA.plot(767,'Beamformed image MLA');
ax(1) = gca;

% beamforming using the "unified pixelbased beamforming" model from 
% Nguyen, N. Q., & Prager, R. W. (2016). High-Resolution Ultrasound Imaging 
% With Unified Pixel-Based Beamforming. IEEE Trans. Med. Imaging, 35(1), 98?108.
mid_MLA_unified_fix =midprocess.das();
mid_MLA_unified_fix.dimension = dimension.both();

mid_MLA_unified_fix.channel_data=channel_data;
mid_MLA_unified_fix.scan=scan_MLA;
mid_MLA_unified_fix.use_unified_fix = 1; %When this flax is set this model is used

mid_MLA_unified_fix.transmit_apodization.window=uff.window.scanline; 
mid_MLA_unified_fix.transmit_apodization.MLA = MLA;
mid_MLA_unified_fix.transmit_apodization.MLA_overlap = MLA/2;

mid_MLA_unified_fix.receive_apodization.window=uff.window.boxcar;
mid_MLA_unified_fix.receive_apodization.f_number=1.7;
b_data_MLA_unified_fix=mid_MLA_unified_fix.go();

b_data_MLA_unified_fix.plot(777,'Beamformed image MLA with unified fix');
ax(2) = gca;

% beamforming using a simpler model assuming PW around focus
mid_MLA_with_plane_fix=midprocess.das();
mid_MLA_with_plane_fix.dimension = dimension.both();
mid_MLA_with_plane_fix.use_PW_fix = 1; % Set this flag to use this model
mid_MLA_with_plane_fix.margin_in_m = 3/1000; %Optionally set the margin of the region around focus to use PW tx delay

mid_MLA_with_plane_fix.channel_data=channel_data;
mid_MLA_with_plane_fix.scan=scan_MLA;

mid_MLA_with_plane_fix.transmit_apodization.window=uff.window.scanline; 
mid_MLA_with_plane_fix.transmit_apodization.MLA = MLA;
mid_MLA_with_plane_fix.transmit_apodization.MLA_overlap = MLA/2;

mid_MLA_with_plane_fix.receive_apodization.window=uff.window.boxcar;
mid_MLA_with_plane_fix.receive_apodization.f_number=1.7;
b_data_MLA_with_plane_fix=mid_MLA_with_plane_fix.go();

b_data_MLA_with_plane_fix.plot(778,'Beamformed image MLA with PW fix');
ax(3) = gca;
linkaxes(ax);

%% Create plot to be used in abstract showing the images and the TX delays
% The images can be zoomed in on the aartifact as we did in the abstract
% We are plotting the TX delay used for the center transmit beam
tx_delay_virtual_source = reshape(mid_MLA.tx_delay_hack,scan_MLA.N_z_axis,scan_MLA.N_x_axis,channel_data.N_waves);
tx_delay_unified_fix = reshape(mid_MLA_unified_fix.tx_delay_hack,scan_MLA.N_z_axis,scan_MLA.N_x_axis,channel_data.N_waves);
tx_delay_plane_fix = reshape(mid_MLA_with_plane_fix.tx_delay_hack,scan_MLA.N_z_axis,scan_MLA.N_x_axis,channel_data.N_waves);

figure(1); clf;
b_data_MLA.plot(subplot(2,3,1),'1a : Virtual source');
ax(1) = gca;
b_data_MLA_unified_fix.plot(subplot(2,3,2),'1b : Model from [2]');
ax(2) = gca;
b_data_MLA_with_plane_fix.plot(subplot(2,3,3),'1c : Virtual source + PW model');
ax(3) = gca;
linkaxes(ax);
subplot(2,3,4); imagesc(scan_MLA.x_axis*1000, scan_MLA.z_axis*1000, tx_delay_virtual_source(:,:,channel_data.N_waves/2));
title('1d: Tx delay virtual source');xlabel('x [mm]');ylabel('z [mm]');colorbar; set(gca,'fontsize',14); 
subplot(2,3,5); imagesc(scan_MLA.x_axis*1000, scan_MLA.z_axis*1000, tx_delay_unified_fix(:,:,channel_data.N_waves/2)); 
title('1e: Tx delay model from [2]');xlabel('x [mm]');ylabel('z [mm]');colorbar; set(gca,'fontsize',14);
subplot(2,3,6); imagesc(scan_MLA.x_axis*1000, scan_MLA.z_axis*1000, tx_delay_plane_fix(:,:,channel_data.N_waves/2));
title('1f: Tx delay virtual source + PW model');xlabel('x [mm]');ylabel('z [mm]');colorbar; set(gca,'fontsize',14);