clear all;
close all;

%% Download data
url='https://nyhirse.medisin.ntnu.no/ustb/data/gcnr/';   % if not found data will be downloaded from here
filename='insilico_20.uff';
tools.download(filename, url, data_path);   

%% Load data
mix = uff.channel_data();
mix.read([data_path filesep filename],'/mix');
channel_SNR = h5read([data_path filesep filename],'/channel_SNR');

%% Scan
sca=uff.linear_scan('x_axis',linspace(-6e-3,6e-3,256).','z_axis', linspace(14e-3,26e-3,256).');

%% Regions

% cyst geometry -> this should go in the uff
x0=0e-3;                
z0=20e-3; 
r=3e-3;                 

% stand off distance <- based on aperture size
M = 55;                             % aperture size
aperture = M * mix.probe.pitch;     % active aperture
F = z0 / aperture;                  % F-number
r_off = round(1.2 * mix.lambda * F, 5); % stand-off distance (rounded to 0.01 mm) 

% boundaries
ri=r-r_off;
ro=r+r_off;
rO=sqrt(ri^2+ro^2);
Ai=pi*ri^2;
Ao=pi*rO^2-pi*ro^2;
d=sqrt((sca.x-x0).^2+(sca.z-z0).^2);

% masks
mask_i=d<ri;
mask_o=(d>ro)&(d<rO);

%% Prepare beamforming
pipe=pipeline();
pipe.scan=sca;
pipe.channel_data = mix;

pipe.transmit_apodization.window=uff.window.flat;
pipe.transmit_apodization.f_number = 1;
pipe.transmit_apodization.minimum_aperture = M*mix.probe.pitch;
pipe.transmit_apodization.maximum_aperture = M*mix.probe.pitch;

pipe.receive_apodization.window=uff.window.flat;
pipe.receive_apodization.f_number = 1;
pipe.receive_apodization.minimum_aperture = M*mix.probe.pitch;
pipe.receive_apodization.maximum_aperture = M*mix.probe.pitch;

das=midprocess.das();

%% DAS

% beamform
das.dimension = dimension.both;
b_das = pipe.go({ das });
%%
close all;
% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast_NLM(M, channel_SNR, b_das, mask_o, mask_i, 'DAS',5);

f1 = figure(1)
hFig1Axes = findobj('Parent',f1,'Type','axes');
f2 = figure(2)
hFig2Axes = findobj('Parent',f2,'Type','axes');
f3 = figure(3)
hFig3Axes = findobj('Parent',f3,'Type','axes');
f4 = figure(4)
hFig4Axes = findobj('Parent',f4,'Type','axes');

images = b_das.get_image();
figure(888);
subplot(531);
imagesc(sca.x_axis*1000,sca.z_axis*1000,images(:,:,7))
caxis([-60 0]); colormap gray;title(['DAS b-mode SNR=',num2str(10*log10(channel_SNR(7)))]);
axis image;
set(gca,'FontSize', 14);
xlabel(gca,'x [mm]');ylabel(gca,'z [mm]');
%viscircles(gca,[x0*1000,z0*1000],ri*1000,'EdgeColor','r','EnhanceVisibility',0);
%viscircles(gca,[x0*1000,z0*1000],ro*1000,'EdgeColor','b','EnhanceVisibility',0);
%viscircles(gca,[x0*1000,z0*1000],rO*1000,'EdgeColor','b','EnhanceVisibility',0);

hSub2 = subplot(532);
copyobj(get(hFig1Axes,'Children'),hSub2);
legend('Field II simulation','Theoretical DAS','Location','se');
axis tight; grid on;
xlabel('||s||');
ylabel('Probability');
set(gca,'FontSize', 14);
legend('p_i','p_o','OVL');
title(['PDFs for SNR=',num2str(10*log10(channel_SNR(7)))]);

% hSub3 = subplot(543);
% copyobj(get(hFig2Axes,'Children'),hSub3);
% legend('Field II simulation','Theoretical DAS');
% axis tight; grid on;
% xlabel('P_{F}');
% ylabel('P_{D}');
% set(gca,'FontSize', 14);
% xlim([0 1])
% ylim([0 1])
% title(['ROC curve for SNR=',num2str(10*log10(channel_SNR(5)))]);

hSub4 = subplot(533);
copyobj(get(hFig4Axes,'Children'),hSub4);
legend('Field II simulation','Theoretical DAS','Location','se');
axis tight; grid on;
xlabel('SNR');ylabel('GCNR');title('DAS GCNR');
set(gca,'FontSize', 14);

%%

b_das_one_img = uff.beamformed_data(b_das);
b_das_one_img.data = b_das_one_img.data(:,:,:,20);
%%


nlm = postprocess.non_local_means_filtering()
% %nlm.flag = 'rician';
nlm.rs = [30 30 1]
nlm.rc = [30 30 1]
nlm.ps = 4;
nlm.sigma = 80;
nlm.input = b_das;

b_das_nlm = nlm.go();


%%
img_db = b_das_one_img.get_image();
clip_at=0;
dynamic_range=60;
%img=uint16(zeros(scan.N_z_axis,scan.N_x_axis));
img=uint16((2^16-1)*min(max((img_db+dynamic_range)./(clip_at+dynamic_range),0),1));

figure;
subplot(131);
imagesc(sca.x_axis*1000,sca.z_axis*1000,img);colormap gray;
axis image
viscircles(gca,[x0*1000,z0*1000],ri*1000,'EdgeColor','r','EnhanceVisibility',0);
viscircles(gca,[x0*1000,z0*1000],ro*1000,'EdgeColor','b','EnhanceVisibility',0);
viscircles(gca,[x0*1000,z0*1000],rO*1000,'EdgeColor','b','EnhanceVisibility',0);


subplot(132);
Kaverage = imfilter(img,fspecial('average',30),'replicate')/255;
imagesc(sca.x_axis*1000,sca.z_axis*1000,Kaverage)
axis image
viscircles(gca,[x0*1000,z0*1000],ri*1000,'EdgeColor','r','EnhanceVisibility',0);
viscircles(gca,[x0*1000,z0*1000],ro*1000,'EdgeColor','b','EnhanceVisibility',0);
viscircles(gca,[x0*1000,z0*1000],rO*1000,'EdgeColor','b','EnhanceVisibility',0);


subplot(133);
imagesc(sca.x_axis*1000,sca.z_axis*1000,b_das_nlm.get_image('none'))
axis image
viscircles(gca,[x0*1000,z0*1000],ri*1000,'EdgeColor','r','EnhanceVisibility',0);
viscircles(gca,[x0*1000,z0*1000],ro*1000,'EdgeColor','b','EnhanceVisibility',0);
viscircles(gca,[x0*1000,z0*1000],rO*1000,'EdgeColor','b','EnhanceVisibility',0);

%%
img_db = b_das.get_image();
clip_at=0;
dynamic_range=60;
for i = 1:20
    i
    %img=uint16((2^16-1)*min(max((img_db(:,:,i)+dynamic_range)./(clip_at+dynamic_range),0),1));
    img_filtered(:,:,i) = double(imfilter(img_db(:,:,i),fspecial('average',30),'replicate'));
    img_filtered(:,:,i) = img_filtered(:,:,i)-max(max(img_filtered(:,:,i)));
    img_filtered(:,:,i) = 10.^(img_filtered(:,:,i)/20);
end

b_mode_average_filt = uff.beamformed_data(b_das); % ToDo: instead we should copy everything but the data

b_mode_average_filt.data = reshape(img_filtered,sca.N_x_axis*sca.N_z_axis,1,1,b_das.N_frames);

[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_mode_average_filt, mask_o, mask_i, 'Average',5);

%%
b_das_nlm.plot
viscircles(gca,[x0*1000,z0*1000],ri*1000,'EdgeColor','r','EnhanceVisibility',0);
viscircles(gca,[x0*1000,z0*1000],ro*1000,'EdgeColor','b','EnhanceVisibility',0);
viscircles(gca,[x0*1000,z0*1000],rO*1000,'EdgeColor','b','EnhanceVisibility',0);

%%
[C, CNR, Pmax, GCNR]=contrast_NLM(M, channel_SNR, b_das_nlm, mask_o, mask_i, 'NLM',5);

f5 = figure(5)
hFig1Axes = findobj('Parent',f5,'Type','axes');
f6 = figure(6)
hFig2Axes = findobj('Parent',f6,'Type','axes');
f7 = figure(7)
hFig3Axes = findobj('Parent',f7,'Type','axes');
f8 = figure(8)
hFig4Axes = findobj('Parent',f8,'Type','axes');

images = b_das_nlm.get_image();
figure(888);
subplot(5,3,4);
imagesc(sca.x_axis*1000,sca.z_axis*1000,images(:,:,7)-max(max(images(:,:,7))))
caxis([-60 0]); colormap gray;title(['NLM b-mode SNR=',num2str(10*log10(channel_SNR(7)))]);
axis image; xlabel('x [mm]'); ylabel('z [mm]')
set(gca,'FontSize', 14);

hSub2 = subplot(5,3,5);
copyobj(get(hFig1Axes,'Children'),hSub2);
axis tight; grid on;
xlabel('||s||');
ylabel('Probability');
set(gca,'FontSize', 14);
legend('p_i','p_o','OVL');
title(['NLM PDFs for SNR=',num2str(10*log10(channel_SNR(7)))]);

% hSub3 = subplot(5,4,7);
% copyobj(get(hFig2Axes,'Children'),hSub3);
% axis tight; grid on;
% xlabel('P_{F}');
% ylabel('P_{D}');
% set(gca,'FontSize', 14);
% axis 'tight';
% legend 'off';
% %xlim([0 1])
% %ylim([0 1])
% title(['NLM ROC curve for SNR=',num2str(10*log10(channel_SNR(5)))]);

hSub4 = subplot(5,3,6);
copyobj(get(hFig4Axes,'Children'),hSub4);
legend('Field II simulation','Theoretical DAS','Location','se');
axis tight; grid on;
xlabel('SNR');ylabel('GCNR');title('NLM GCNR');
set(gca,'FontSize', 14);


%% S-DAS

% beamform
das.dimension = dimension.both;

slope=10;
xA=[linspace(0,-25,3) -60+linspace(0,25,3)].';
yA=[linspace(0,-25,3)/slope -60+linspace(0,25,3)/slope].';

f=fit(xA,yA,'smoothingspline')
figure;
xx=linspace(-60,0,100); 
plot(xx,xx,'r--','linewidth',2); hold on; grid on; axis equal;
plot(xx,f(xx)','linewidth',2);
xlim([-60 0]);
ylim([-62 2]);
set(gca,'FontSize', 12);
xlabel('Input [dB]');
ylabel('Output [sB]');
legend('linear','S-curve','Location','NorthWest')

% anechoic
b_sdas = uff.beamformed_data(b_das);
b_sdas.data = 20*log10(abs(b_das.data)./max(abs(b_das.data)));
b_sdas.data = 10.^(reshape(f(b_sdas.data),size(b_sdas.data))/20);

% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_sdas, mask_o, mask_i, 'S-DAS');

%% CF

% beamform
das.dimension = dimension.transmit;
cf = postprocess.coherence_factor();
cf.dimension = dimension.receive;
cf_anechoic = pipe.go({ das cf });

%%
% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, cf_anechoic, mask_o, mask_i, 'CF');

f5 = figure(5)
hFig1Axes = findobj('Parent',f5,'Type','axes');
f6 = figure(6)
hFig2Axes = findobj('Parent',f6,'Type','axes');
f7 = figure(7)
hFig3Axes = findobj('Parent',f7,'Type','axes');
f8 = figure(8)
hFig4Axes = findobj('Parent',f8,'Type','axes');

images = cf_anechoic.get_image();
figure(888);
subplot(535);
imagesc(sca.x_axis*1000,sca.z_axis*1000,images(:,:,20)-max(max(images(:,:,20))))
caxis([-60 0]); colormap gray;title(['CF b-mode SNR=',num2str(10*log10(channel_SNR(20)))]);
axis image; axis image; xlabel('x [mm]'); ylabel('z [mm]')
set(gca,'FontSize', 14);

hSub2 = subplot(546);
copyobj(get(hFig1Axes,'Children'),hSub2);
axis tight; grid on;
xlabel('||s||');
ylabel('Probability');
set(gca,'FontSize', 14);
legend('p_i','p_o','OVL');
title(['CF PDFs for SNR=',num2str(10*log10(channel_SNR(5)))]);

hSub3 = subplot(547);
copyobj(get(hFig2Axes,'Children'),hSub3);
axis tight; grid on;
xlabel('P_{F}');
ylabel('P_{D}');
set(gca,'FontSize', 14);
axis 'tight';
legend 'off';
%xlim([0 1])
%ylim([0 1])
title(['CF ROC curve for SNR=',num2str(10*log10(channel_SNR(5)))]);

hSub4 = subplot(548);
copyobj(get(hFig4Axes,'Children'),hSub4);
legend('Field II simulation','Theoretical DAS','Location','se');
axis tight; grid on;
xlabel('SNR');ylabel('GCNR');title('CF GCNR');
set(gca,'FontSize', 14);



%%
nlm = postprocess.non_local_means_filtering();
nlm.input = cf_anechoic;

b_cf_nlm = nlm.go();

[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_cf_nlm, mask_o, mask_i, 'NLM');
%% PCF

% beamform
das.dimension = dimension.transmit;
pcf = postprocess.phase_coherence_factor();
pcf.center_frequency = 5e6;
pcf.dimension = dimension.receive;
pcf.gamma=1;
pipe.channel_data=mix;
pcf_anechoic = pipe.go({ das pcf });

%% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast_NLM(M, channel_SNR, pcf_anechoic, mask_o, mask_i, 'PCF');


f5 = figure(9)
hFig1Axes = findobj('Parent',f5,'Type','axes');
f6 = figure(10)
hFig2Axes = findobj('Parent',f6,'Type','axes');
f7 = figure(11)
hFig3Axes = findobj('Parent',f7,'Type','axes');
f8 = figure(12)
hFig4Axes = findobj('Parent',f8,'Type','axes');

images = pcf_anechoic.get_image();
figure(888);
subplot(5,3,7);
imagesc(sca.x_axis*1000,sca.z_axis*1000,images(:,:,5)-max(max(images(:,:,5))))
caxis([-60 0]); colormap gray;title(['PCF b-mode SNR=',num2str(10*log10(channel_SNR(5)))]);
axis image; xlabel('x [mm]'); ylabel('z [mm]')
set(gca,'FontSize', 14);

hSub2 = subplot(5,3,8);
copyobj(get(hFig1Axes,'Children'),hSub2);
axis tight; grid on;
xlabel('||s||');
ylabel('Probability');
set(gca,'FontSize', 14);
legend('p_i','p_o','OVL');
title(['PCF PDFs for SNR=',num2str(10*log10(channel_SNR(7)))]);

% hSub3 = subplot(5,3,9);
% copyobj(get(hFig2Axes,'Children'),hSub3);
% axis tight; grid on;
% xlabel('P_{F}');
% ylabel('P_{D}');
% set(gca,'FontSize', 14);
% axis 'tight';
% legend 'off';
% %xlim([0 1])
% %ylim([0 1])
% title(['PCF ROC curve for SNR=',num2str(10*log10(channel_SNR(5)))]);

hSub4 = subplot(5,3,9);
copyobj(get(hFig4Axes,'Children'),hSub4);
legend('Field II simulation','Theoretical DAS','Location','se');
axis tight; grid on;
xlabel('SNR');ylabel('GCNR');title('PCF GCNR');
set(gca,'FontSize', 14);


%% GCF

% beamform
das.dimension = dimension.transmit;
gcf = postprocess.generalized_coherence_factor;
gcf.dimension = dimension.receive;
gcf.M0=4;
pipe.channel_data=mix;
gcf_anechoic = pipe.go({ das gcf });
%%

% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, gcf.GCF, mask_o, mask_i, 'GCF');


f5 = figure(13)
hFig1Axes = findobj('Parent',f5,'Type','axes');
f6 = figure(14)
hFig2Axes = findobj('Parent',f6,'Type','axes');
f7 = figure(15)
hFig3Axes = findobj('Parent',f7,'Type','axes');
f8 = figure(16)
hFig4Axes = findobj('Parent',f8,'Type','axes');

images = gcf_anechoic.get_image();
figure(888);
subplot(5,3,10);
imagesc(sca.x_axis*1000,sca.z_axis*1000,images(:,:,20)-max(max(images(:,:,20))))
caxis([-60 0]); colormap gray;title(['GCF b-mode SNR=',num2str(10*log10(channel_SNR(20)))]);
axis image; xlabel('x [mm]'); ylabel('z [mm]')
set(gca,'FontSize', 14);

hSub2 = subplot(5,3,11);
copyobj(get(hFig1Axes,'Children'),hSub2);
axis tight; grid on;
xlabel('||s||');
ylabel('Probability');
set(gca,'FontSize', 14);
legend('p_i','p_o','OVL');
title(['GCF PDFs for SNR=',num2str(10*log10(channel_SNR(5)))]);

% hSub3 = subplot(5,4,15);
% copyobj(get(hFig2Axes,'Children'),hSub3);
% axis tight; grid on;
% xlabel('P_{F}');
% ylabel('P_{D}');
% set(gca,'FontSize', 14);
% axis 'tight';
% legend 'off';
% %xlim([0 1])
% %ylim([0 1])
% title(['GCF ROC curve for SNR=',num2str(10*log10(channel_SNR(5)))]);

hSub4 = subplot(5,3,12);
copyobj(get(hFig4Axes,'Children'),hSub4);
legend('Field II simulation','Theoretical DAS','Location','se');
axis tight; grid on;
xlabel('SNR');ylabel('GCNR');title('GCF GCNR');
set(gca,'FontSize', 14);

%%

%% DMAS

% % beamform
% das.dimension = dimension.transmit;
% b_data_tx = pipe.go({das});
% 
% % process DMAS
% dmas=postprocess.delay_multiply_and_sum();
% dmas.dimension = dimension.receive;
% dmas.transmit_apodization = pipe.transmit_apodization;
% dmas.receive_apodization = pipe.receive_apodization;
% dmas.input = b_data_tx;
% dmas.channel_data = mix;
% b_dmas = dmas.go();
% 
% % evaluate contrast
% [C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_dmas, mask_o, mask_i, 'DMAS');


%% SLSC using M elements (hacking the SLSC postprocess)

% important that we use only M elements, centered around the abscissa of the pixel. 
% Changing that will alter the SNR ratio.
das.dimension = dimension.transmit;
b_data_tx = pipe.go({das});

%%
% Set up the SLSC postprocess
slsc = postprocess.short_lag_spatial_coherence();
slsc.receive_apodization = das.receive_apodization;
slsc.dimension = dimension.receive;
slsc.channel_data = mix;
slsc.maxM = 14;
slsc.input = b_data_tx;
slsc.K_in_lambda = 1;

Q = slsc.maxM./M

% Pick out the M channels that contains data, thus only the M elemetns
% centered around the abscissa of the pixel. This is taken care of by the
% apodization so we can just pick up the element signals that are larger than 0.

% Data buffers
aux_data = zeros(b_data_tx.scan.N_z_axis*b_data_tx.scan.N_x_axis,1,1,mix.N_frames);
aux_data_clamped = zeros(b_data_tx.scan.N_z_axis*b_data_tx.scan.N_x_axis,1,1,mix.N_frames);
data_cube_M_elements = complex(zeros(b_data_tx.scan.N_z_axis,b_data_tx.scan.N_x_axis,M,1));
for f = 1:mix.N_frames
    f
    % Reshape the beamformed data as a cube (Z,X,Elements)
    data_cube = reshape(b_data_tx.data(:,:,1,f),sca.N_z_axis,sca.N_x_axis,mix.N_channels);
    for x = 1: b_data_tx.scan.N_x_axis
        sum_over_z = abs(sum(squeeze(data_cube(:,x,:)),1));
        elements_with_data = sum_over_z>0;
        data_cube_M_elements(:,x,:) = data_cube(:,x,elements_with_data);
    end
    
    % Call the actual implementation of the SLSC calculations, but using
    % the cube which has only M elements
    [image,slsc_values] = slsc.short_lag_spatial_coherence_implementation(data_cube_M_elements);
    image(image<0) = 0; % Set negative coherence values to zero
    
    slsc_img = squeeze(sum(slsc_values(:,:,:),2));
    slsc_img = slsc_img./max(slsc_img(:)); %According to previous publications
    
    % Make one clamped version
    slsc_img_clamped = slsc_img;
    slsc_img_clamped(slsc_img_clamped < 0 ) = 0;
    aux_data_clamped(:,1,1,f) = slsc_img_clamped(:);
    
    % In the other we  can shift the negative coherence to something
    % positive.
    slsc_img = slsc_img + abs(min(slsc_img(:)));
    aux_data(:,1,1,f) = slsc_img(:);
end

% Put the resulting SLSC images in a beamformed data
b_slsc_M_shifted = uff.beamformed_data();
b_slsc_M_shifted.scan = sca;
b_slsc_M_shifted.data = aux_data;

b_slsc_M_clamped = uff.beamformed_data();
b_slsc_M_clamped.scan = sca;
b_slsc_M_clamped.data = aux_data_clamped;

%%
%[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_slsc_M_shifted, mask_o, mask_i, 'SLSC');
idx = 1;
    close all;
    to_plot = 5;
    text_size = 13;
for colms = 1:6

    
    if colms == 1
        tag = 'DAS';
        [C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_das, mask_o, mask_i, tag, to_plot);
        images = b_das.get_image();
        f5 = figure(1)
        hFig1Axes = findobj('Parent',f5,'Type','axes');
        f6 = figure(2)
        hFig2Axes = findobj('Parent',f6,'Type','axes');
        f7 = figure(3)
        hFig3Axes = findobj('Parent',f7,'Type','axes');
        f8 = figure(4)
        hFig4Axes = findobj('Parent',f8,'Type','axes');
        
    elseif colms == 2
        tag = 'NLM';
        [C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_das_nlm, mask_o, mask_i, tag, to_plot);
        images = b_das_nlm.get_image();
        
        f5 = figure(5)
        hFig1Axes = findobj('Parent',f5,'Type','axes');
        f6 = figure(6)
        hFig2Axes = findobj('Parent',f6,'Type','axes');
        f7 = figure(7)
        hFig3Axes = findobj('Parent',f7,'Type','axes');
        f8 = figure(8)
        hFig4Axes = findobj('Parent',f8,'Type','axes');
        
    elseif colms == 3
        %%
        tag = 'GF'
        [C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_mode_average_filt, mask_o, mask_i, tag,to_plot);
        images = b_mode_average_filt.get_image()
        
        f5 = figure(9)
        hFig1Axes = findobj('Parent',f5,'Type','axes');
        f6 = figure(10)
        hFig2Axes = findobj('Parent',f6,'Type','axes');
        f7 = figure(11)
        hFig3Axes = findobj('Parent',f7,'Type','axes');
        f8 = figure(12)
        hFig4Axes = findobj('Parent',f8,'Type','axes');
    elseif colms == 4
        tag = 'PCF';
        [C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, pcf_anechoic, mask_o, mask_i, tag, to_plot);
        images = pcf_anechoic.get_image();
        
        f5 = figure(14)
        hFig1Axes = findobj('Parent',f5,'Type','axes');
        f6 = figure(14)
        hFig2Axes = findobj('Parent',f6,'Type','axes');
        f7 = figure(15)
        hFig3Axes = findobj('Parent',f7,'Type','axes');
        f8 = figure(16)
        hFig4Axes = findobj('Parent',f8,'Type','axes');
        
    elseif colms == 5
        tag = 'GCF';
        [C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, gcf_anechoic, mask_o, mask_i, tag, to_plot);
        images = gcf_anechoic.get_image();
        
        f5 = figure(17)
        hFig1Axes = findobj('Parent',f5,'Type','axes');
        f6 = figure(18)
        hFig2Axes = findobj('Parent',f6,'Type','axes');
        f7 = figure(19)
        hFig3Axes = findobj('Parent',f7,'Type','axes');
        f8 = figure(20)
        hFig4Axes = findobj('Parent',f8,'Type','axes');
        
    elseif colms == 6
        tag = 'SLSC';
        [C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_slsc_M_clamped, mask_o, mask_i, tag, to_plot);
        images = b_slsc_M_clamped.get_image();
        
        f5 = figure(21)
        hFig1Axes = findobj('Parent',f5,'Type','axes');
        f6 = figure(22)
        hFig2Axes = findobj('Parent',f6,'Type','axes');
        f7 = figure(23)
        hFig3Axes = findobj('Parent',f7,'Type','axes');
        f8 = figure(24)
        hFig4Axes = findobj('Parent',f8,'Type','axes');
        
    end

    
    
    figure(888);
    subplot(3,6,colms);
    imagesc(sca.x_axis*1000,sca.z_axis*1000,images(:,:,20)-max(max(images(:,:,20))))
    caxis([-60 0]); colormap gray;title([tag,' SNR=',num2str(10*log10(channel_SNR(20)),'%.1f')]);
    axis image; xlabel('x [mm]'); ylabel('z [mm]')
    set(gca,'FontSize', text_size);
    idx = idx + 1;
    colorbar;
    
    hSub2 = subplot(3,6,colms+6);
    copyobj(get(hFig2Axes,'Children'),hSub2);
    axis tight; grid on;
    xlabel('||s||');
    ylabel('Probability');
    set(gca,'FontSize', text_size);
    if colms == 1
        legend('p_i','p_o','OVL');
    end
    title([tag,' SNR=',num2str(10*log10(channel_SNR(to_plot)),'%.1f')]);
    idx = idx + 1;
    if colms == 4
        %xlim([0 10])
        xlim([0 5])
    elseif colms == 5
        %xlim([0 20])
        xlim([0 30])
    elseif colms == 6
        %xlim([0 0.2])
        %xlim([0 0.3])
        ylim([0 0.02])
    end
    
    hSub4 = subplot(3,6,colms+12);
    copyobj(get(hFig4Axes,'Children'),hSub4);
    if colms == 1
        legend('Field II simulation','Theoretical DAS','Location','se');
    end
    axis tight; grid on;
    xlabel('SNR');ylabel('GCNR');title([tag,' GCNR']);
    set(gca,'FontSize', text_size);
    idx = idx + 1;
end
%%

%%

nlm = postprocess.non_local_means_filtering();
nlm.input = b_slsc_M_clamped;
nlm.run_on_logcompressed = 0;
nlm.sigma = 100;
nlm.beta = 1.0;     % default 1.0
nlm.rs   = [10,10,1]; % default [2,2,2]
nlm.rc   = [5,5,1]; % default [1,1,1]
nlm.ps   = 0.5;       % default 2
nlm.flag = 'gaussian'; % default 'gaussian', could be rician
nlm.block= 1;       % 0 loop, 1:vector operations (for larger search windows)


b_SLSC_NLM = nlm.go();
%%
% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_SLSC_NLM, mask_o, mask_i, 'NLM on SLSC');


f5 = figure(13)
hFigIAxes = findobj('Parent',f5,'Type','axes');
f6 = figure(14)
hFig2Axes = findobj('Parent',f6,'Type','axes');
f7 = figure(15)
hFig3Axes = findobj('Parent',f7,'Type','axes');
f8 = figure(16)
hFig4Axes = findobj('Parent',f8,'Type','axes');

images = b_SLSC_NLM.get_image();
figure(888);
subplot(4,4,13);
imagesc(sca.x_axis*1000,sca.z_axis*1000,images(:,:,20))
caxis([-60 0]); colormap gray;title(['NLM on SLSC SNR=',num2str(10*log10(channel_SNR(end)))]);
axis image;
hSub2 = subplot(4,4,14);
copyobj(get(hFig2Axes,'Children'),hSub2);
legend('Field II simulation','Theoretical DAS');
axis tight;grid on
xlabel('SNR');ylabel('C');title('NLM on SLSC C');

hSub3 = subplot(4,4,15);
copyobj(get(hFig3Axes,'Children'),hSub3);
legend('Field II simulation','Theoretical DAS','Location','se');
axis tight;grid on
xlabel('SNR');ylabel('CNR');title('NLM on SLSC CNR');

hSub4 = subplot(4,4,16);
copyobj(get(hFig4Axes,'Children'),hSub4);
legend('Field II simulation','Theoretical DAS','Location','se');
axis tight;grid on
xlabel('SNR');ylabel('GCNR');title('NLM on SLSC GCNR');



%% SLSC with full aperture
das.dimension = dimension.transmit;
pipe.transmit_apodization.window = uff.window.none;
pipe.receive_apodization.window = uff.window.none;

% Set up the SLSC postprocess
slsc = postprocess.short_lag_spatial_coherence();
slsc.receive_apodization = pipe.receive_apodization;
slsc.dimension = dimension.receive;
slsc.channel_data = mix;
slsc.maxM = 20;
slsc.K_in_lambda = 1;

b_slsc = pipe.go({das slsc});

%% evaluate contrast
[C, CNR, Pmax, GCNR]=contrast(M, channel_SNR, b_slsc, mask_o, mask_i, 'SLSC (N=20) full aperture');


