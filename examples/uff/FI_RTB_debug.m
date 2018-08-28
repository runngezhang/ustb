%% Focused Linear scan with L7-4 probe in Verasonics demonstrating RTB
%
% This script is available in the USTB repository as
% examples/uff/FI_UFF_Verasonics_linear_RTB.m
%
% This example demonstrates the RTB implementation and demonstrates
% different fixes to the artifact occuring near the focus.
%
% One solution is the transmit delay model introduced in  Nguyen, N. Q., &
% Prager, R. W. (2016). High-Resolution Ultrasound Imaging With Unified Pixel-Based 
% Beamforming. IEEE Trans. Med. Imaging, 35(1), 98-108.
%
% Another solution is a simpler model assuming PW around focus.
%
% This scripts creates the figure used in the abstract submitted to
% the IEEE IUS 2018 with title "A simple, artifact-free virtual source model"
%   
% _by Ole Marius Hoel Rindal <olemarius@olemarius.net> Last updated: 2018/10/05

%% Read channel data

clear all; close all;

% data location
url='http://ustb.no/datasets/';      % if not found downloaded from here

filename='L7_FI_IUS2018.uff';
tools.download(filename, url, data_path);   
channel_data = uff.read_object([data_path filesep filename],'/channel_data');

%%
    figure(1);
    
for i = 1:128
    subplot(211);clf;hold on;
    plot(channel_data.data(:,1,1),'b');
    plot(channel_data.data(:,i,1),'r');
    drawnow;
    pause(0.5);
end

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

%%
% Display image
b_data.plot([],'Conventional one scanline per transmit');


%% Retrospective beamforming (RTB) with conventional virtual source model
% Create scan with MLA's
MLA = 4;
scan_RTB=uff.linear_scan('x_axis',linspace(x_axis(1),x_axis(end),...
                                    length(x_axis)*MLA)','z_axis',z_axis);

%% RTB using a simpler model assuming PW around focus
mid_RTB_with_plane_fix=midprocess.das();
mid_RTB_with_plane_fix.dimension = dimension.receive();
mid_RTB_with_plane_fix.transmit_delay_model = transmit_delay_model.hybrid;
%Optionally set the margin of the region around focus to use PW tx delay
mid_RTB_with_plane_fix.pw_margin = 1/1000; 

mid_RTB_with_plane_fix.channel_data=channel_data;
mid_RTB_with_plane_fix.scan=scan_RTB;

mid_RTB_with_plane_fix.transmit_apodization.window=uff.window.tukey50;
mid_RTB_with_plane_fix.transmit_apodization.f_number = 3;
mid_RTB_with_plane_fix.transmit_apodization.MLA = MLA;
mid_RTB_with_plane_fix.transmit_apodization.MLA_overlap = 1;

mid_RTB_with_plane_fix.receive_apodization.window=uff.window.boxcar;
mid_RTB_with_plane_fix.receive_apodization.f_number=1.7;
b_data_RTB_with_plane_fix=mid_RTB_with_plane_fix.go();

% Calculate the transmit apodzation used to compensate image
tx_apod = mid_RTB_with_plane_fix.transmit_apodization.data;

b_data_RTB_plane_fix_weighted = uff.beamformed_data(b_data_RTB_with_plane_fix);
b_data_RTB_plane_fix_weighted.data = b_data_RTB_plane_fix_weighted.data...
                                                        .*(1./sum(tx_apod,2));
b_data_RTB_plane_fix_weighted.plot(10,'RTB image with PW hybrid virtual source model');

%%
images = reshape(b_data_RTB_plane_fix_weighted.data,[scan_RTB.N_z_axis scan_RTB.N_x_axis 128]);


figure(888);
subplot(221);
imagesc(db(abs(images(:,:,1))));colorbar;caxis([-50 80])
subplot(222);
imagesc(db(abs(images(:,:,128))));colorbar;caxis([-50 80])


%%

img_incoherent = sum(abs(images),3);
img_coherent = sum(images,3);
figure(44);
subplot(131)
imagesc(b_data.get_image);caxis([-60 0]);colormap gray;title('scanline');
subplot(133)
imagesc(db(img_coherent./max(img_coherent(:))));caxis([-60 0]);colormap gray;title('coherent')
subplot(132)
imagesc(db(img_incoherent./max(img_incoherent(:))));caxis([-60 0]);colormap gray;title('incoherent');
%%
figure(88);
subplot(211)
imagesc(db(abs(images(:,:,116))))
subplot(212)
imagesc(db(abs(images(:,:,117))))

figure(89);
subplot(211)
imagesc(db(abs(images(:,:,21))))
subplot(212)
imagesc(db(abs(images(:,:,22))))
%%

figure(89);clf;
subplot(221);hold on;
plot(real(images(:,467,116)),'b','DisplayName','transmit 116');ylim([-1000 1000])
plot(real(images(:,467,117)),'r','DisplayName','transmit 117');ylim([-1000 1000])
title('Real part of axial line 467');
ax(1) = gca;
subplot(222);hold on;
plot(imag(images(:,467,116)),'b','DisplayName','transmit 116');ylim([-1000 1000])
plot(imag(images(:,467,117)),'r','DisplayName','transmit 117');ylim([-1000 1000])
title('Imag part of axial line 467');
ax(2) = gca;
subplot(223);hold on;
plot(real(images(:,80,21)),'b','DisplayName','transmit 21');ylim([-2000 2000])
plot(real(images(:,80,22)),'r','DisplayName','transmit 22');ylim([-2000 2000])
title('Real part of axial line 80');
ax(3) = gca;
subplot(224);hold on;
plot(imag(images(:,80,21)),'b','DisplayName','transmit 21');ylim([-2000 2000])
plot(imag(images(:,80,22)),'r','DisplayName','transmit 22');ylim([-2000 2000])
title('Imag part of axial line 80');
ax(4) = gca;
linkaxes(ax);

%%
figure(91);clf;
subplot(121);hold on;
plot(abs(images(:,467,116)),'b');ylim([-1000 1000])
plot(abs(images(:,467,117)),'r');ylim([-1000 1000])
ax(1) = gca;
subplot(122);hold on;
plot(abs(images(:,80,21)),'b');
plot(abs(images(:,80,22)),'r');
ax(2) = gca;
linkaxes(ax);

%%
figure(66);
subplot(121);hold on;
plot(channel_data.data(:,100,116),'b');ylim([-10000 10000])
plot(channel_data.data(:,100,117),'r');
ax(1) = gca;
subplot(122);hold on;
plot(channel_data.data(:,20,21),'b');
plot(channel_data.data(:,20,22),'r');
ax(2) = gca;
linkaxes(ax);