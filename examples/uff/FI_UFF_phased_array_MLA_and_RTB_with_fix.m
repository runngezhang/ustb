% Demonstrating MLA and RTB for a phased array

% Clear up
clear all;close all;

% Read the data, poentitally download it
url='http://ustb.no/datasets/';      % if not found downloaded from here
local_path = [ustb_path(),'/data/']; % location of example data
addpath(local_path);

% Choose dataset
filename='P4_FI.uff';
% check if the file is available in the local path or downloads otherwise
tools.download(filename, url, local_path);
channel_data = uff.read_object([local_path, filename],'/channel_data');

%% Throw away most of the dataset to make the beamformign run faster. 
% Comment out this line if you want to run it on all 50 frames
channel_data.N_frames = 1;

%% Create the sector scan we want to reconstruct
scan=uff.sector_scan('azimuth_axis',...
    linspace(channel_data.sequence(1).source.azimuth,channel_data.sequence(end).source.azimuth,...
    length(channel_data.sequence))'...
    ,'depth_axis',linspace(0,90e-3,512)');

%%
mid = midprocess.das();
mid.channel_data=channel_data;
mid.dimension = dimension.both();
mid.scan=scan;
mid.transmit_apodization.window=uff.window.scanline;
mid.receive_apodization.window=uff.window.tukey25;
mid.receive_apodization.f_number = 1.7;

b_data = mid.go();
%% Plot the image 
b_data.plot(40,['DAS']);
ax(1) = gca;

%% Beamform the image with 4 MLA's per scan line with two overlapping
MLA = 8;
scan_MLA=uff.sector_scan('azimuth_axis',...
    linspace(scan.azimuth_axis(1),scan.azimuth_axis(end),scan.N_azimuth_axis*MLA)'...
    ,'depth_axis',scan.depth_axis);

%%
mid_MLA=midprocess.das();
mid_MLA.channel_data=channel_data;
mid_MLA.dimension = dimension.both();
mid_MLA.scan=scan_MLA;
mid_MLA.transmit_apodization.window=uff.window.scanline;
mid_MLA.transmit_apodization.MLA = MLA;
mid_MLA.transmit_apodization.MLA_overlap = MLA;
mid_MLA.receive_apodization.window=uff.window.tukey25;
mid_MLA.receive_apodization.f_number = 1.7;

b_data_MLA = mid_MLA.go();
%% Plot the image 
b_data_MLA.plot(41,['DAS with MLAs']);
ax(2) = gca;

%%

mid_MLA_unified_fix=midprocess.das();
mid_MLA_unified_fix.channel_data=channel_data;
mid_MLA_unified_fix.dimension = dimension.both();
mid_MLA_unified_fix.use_unified_fix = 1;
mid_MLA_unified_fix.scan=scan_MLA;
mid_MLA_unified_fix.transmit_apodization.window=uff.window.scanline;
mid_MLA_unified_fix.transmit_apodization.MLA = MLA;
mid_MLA_unified_fix.transmit_apodization.MLA_overlap = MLA;
mid_MLA_unified_fix.receive_apodization.window=uff.window.tukey25;
mid_MLA_unified_fix.receive_apodization.f_number = 1.7;

b_data_MLA_unified_fix = mid_MLA_unified_fix.go();
%% Plot the image 
b_data_MLA_unified_fix.plot(42,['DAS with MLAs unified FIX']);
ax(4) = gca;

linkaxes(ax)
%%
mid_MLA_plane_fix=midprocess.das();
mid_MLA_plane_fix.channel_data=channel_data;
mid_MLA_plane_fix.dimension = dimension.both();
mid_MLA_plane_fix.use_PW_fix = 1;
mid_MLA_plane_fix.margin_in_m = 2;
mid_MLA_plane_fix.scan=scan_MLA;
mid_MLA_plane_fix.transmit_apodization.window=uff.window.scanline;
mid_MLA_plane_fix.transmit_apodization.MLA = MLA;
mid_MLA_plane_fix.transmit_apodization.MLA_overlap = MLA;
mid_MLA_plane_fix.receive_apodization.window=uff.window.tukey25;
mid_MLA_plane_fix.receive_apodization.f_number = 1.7;

b_data_MLA_plane_fix = mid_MLA_plane_fix.go();
%% Plot the image 
b_data_MLA_plane_fix.plot(43,['DAS with MLAs using PW fix']);
ax(5) = gca;

linkaxes(ax)


%% Lets have a closer look at the focal region
% Both regions fix the focal region in front of the tranducer
f100 = figure(100);
set(f100,'Position',[141 284 1256 514]);
b_data.plot(subplot(2,4,1),['DAS']);
ax_sub_top(1) = gca;
b_data_MLA.plot(subplot(2,4,2),['DAS with MLAs']);
ax_sub_top(2) = gca;
b_data_MLA_unified_fix.plot(subplot(2,4,3),['DAS with MLAs unified FIX']);
ax_sub_top(3) = gca;
b_data_MLA_plane_fix.plot(subplot(2,4,4),['DAS with MLAs using PW fix']);
ax_sub_top(4) = gca;
linkaxes(ax_sub_top)
ylim([45 60]);xlim([-10 10]);

% However, there must be something wrong with my implementation of the uniform
% fix to the side
b_data.plot(subplot(2,4,5),['DAS']);
ax_sub_bottom(1) = gca;
b_data_MLA.plot(subplot(2,4,6),['DAS with MLAs']);
ax_sub_bottom(2) = gca;
b_data_MLA_unified_fix.plot(subplot(2,4,7),['DAS with MLAs unified FIX']);
ax_sub_bottom(3) = gca;
b_data_MLA_plane_fix.plot(subplot(2,4,8),['DAS with MLAs using PW fix']);
ax_sub_bottom(4) = gca;
linkaxes(ax_sub_bottom)
ylim([40 55]);xlim([15 35]);

%% Let's try to do full retrospective beamforming
mid_RTB=midprocess.das();
mid_RTB.channel_data=channel_data;
mid_RTB.dimension = dimension.both();
mid_RTB.scan=scan_MLA;
mid_RTB.transmit_apodization.window = uff.window.sector_scan_rtb;
mid_RTB.transmit_apodization.probe = channel_data.probe;
mid_RTB.receive_apodization.window = uff.window.tukey25;
mid_RTB.receive_apodization.f_number = 1.7;

b_data_RTB = mid_RTB.go();

%% Get the transmit apod
tx_apod = mid_RTB.transmit_apodization.data;

x_matrix=reshape(scan_MLA.x,[scan_MLA.N_depth_axis scan_MLA.N_azimuth_axis]);
z_matrix=reshape(scan_MLA.z,[scan_MLA.N_depth_axis scan_MLA.N_azimuth_axis]);

figure(88);
pcolor(x_matrix*1e3,z_matrix*1e3,reshape(tx_apod(:,51),scan_MLA.N_depth_axis,scan_MLA.N_azimuth_axis))
xlabel('x [mm]');
ylabel('z [mm]');
shading('flat');
set(gca,'fontsize',14);
set(gca,'YDir','reverse');
axis('tight','equal');
title('TX apod from sequence 51');

% Calculate weights based on the transmit apod
weighting = 1./sum(tx_apod,2);

b_data_RTB_weighted = uff.beamformed_data(b_data_RTB);
b_data_RTB_weighted.data = b_data_RTB_weighted.data.*weighting;

b_data_RTB_weighted.plot(11,'RTB')
ax(5) = gca;

linkaxes(ax)

%% Let's try to do full retrospective beamforming
mid_RTB_unified_fix=midprocess.das();
mid_RTB_unified_fix.channel_data=channel_data;
mid_RTB_unified_fix.dimension = dimension.both();
mid_RTB_unified_fix.use_unified_fix = 1;
mid_RTB_unified_fix.scan=scan_MLA;
mid_RTB_unified_fix.transmit_apodization.window = uff.window.sector_scan_rtb;
mid_RTB_unified_fix.transmit_apodization.probe = channel_data.probe;
mid_RTB_unified_fix.receive_apodization.window = uff.window.tukey25;
mid_RTB_unified_fix.receive_apodization.f_number = 1.7;

b_data_RTB_unified_fix = mid_RTB_unified_fix.go();

b_data_RTB_unified_fix_weighted = uff.beamformed_data(b_data_RTB_unified_fix);
b_data_RTB_unified_fix_weighted.data = b_data_RTB_unified_fix_weighted.data.*weighting;

b_data_RTB_unified_fix_weighted.plot(12,'RTB with unified fix')
ax(5) = gca;

linkaxes(ax)

%%
mid_RTB_PW_fix=midprocess.das();
mid_RTB_PW_fix.channel_data=channel_data;
mid_RTB_PW_fix.dimension = dimension.both();
mid_RTB_PW_fix.use_PW_fix = 1;
mid_RTB_PW_fix.margin_in_m = 3;
mid_RTB_PW_fix.scan=scan_MLA;
mid_RTB_PW_fix.transmit_apodization.window = uff.window.sector_scan_rtb;
mid_RTB_PW_fix.transmit_apodization.probe = channel_data.probe;
mid_RTB_PW_fix.receive_apodization.window = uff.window.tukey25;
mid_RTB_PW_fix.receive_apodization.f_number = 1.7;

b_data_RTB_PW_fix = mid_RTB_PW_fix.go();

b_data_RTB_PW_fix_weighted = uff.beamformed_data(b_data_RTB_PW_fix);
b_data_RTB_PW_fix_weighted.data = b_data_RTB_PW_fix_weighted.data.*weighting;

b_data_RTB_PW_fix_weighted.plot(13,'RTB with PW fix')
ax(5) = gca;

linkaxes(ax)

%% Lets have a closer look at the focal region
% Both regions fix the focal region in front of the tranducer
f200 = figure(200);
set(f100,'Position',[141 284 1256 514]);
b_data.plot(subplot(2,4,1),['DAS']);
ax_sub_top(1) = gca;
b_data_RTB_weighted.plot(subplot(2,4,2),['RTB']);
ax_sub_top(2) = gca;
b_data_RTB_unified_fix_weighted.plot(subplot(2,4,3),['RTB unified FIX']);
ax_sub_top(3) = gca;
b_data_RTB_PW_fix_weighted.plot(subplot(2,4,4),['RTB using PW fix']);
ax_sub_top(4) = gca;
linkaxes(ax_sub_top)
ylim([45 60]);xlim([-10 10]);

% However, there must be something wrong with my implementation of the uniform
% fix to the side
b_data.plot(subplot(2,4,5),['DAS']);
ax_sub_bottom(1) = gca;
b_data_RTB_weighted.plot(subplot(2,4,6),['RTB']);
ax_sub_bottom(2) = gca;
b_data_RTB_unified_fix_weighted.plot(subplot(2,4,7),['RTB unified FIX']);
ax_sub_bottom(3) = gca;
b_data_RTB_PW_fix_weighted.plot(subplot(2,4,8),['RTB using PW fix']);
ax_sub_bottom(4) = gca;
linkaxes(ax_sub_bottom)
ylim([40 55]);xlim([15 35]);