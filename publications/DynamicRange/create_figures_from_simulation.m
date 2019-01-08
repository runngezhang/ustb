%% Create the figures from the simulated dynamic range phantom
% For the publication Rindal, O. M. H., Austeng, A., Fatemi, A., 
% & Rodriguez-Molares, A. (2018). The effect of dynamic range transformations
% in the estimation of contrast. Submitted to IEEE Transactions on Ultrasonics,
% Ferroelectrics, and Frequency Control.
%
% NB! You need to run the create_figures_from_experimental.m first since we
% are using the stored experimental gradient values and contrast values to
% plot both the simulated and the experimental in the same plot.
%
% Author: Ole Marius Hoel Rindal <olemarius@olemarius.net> 05.06.18

%% Load the data from the uff file
clear all;
close all;

filename = [data_path,filesep,'FieldII_STAI_dynamic_range_similar_to_exp_v5_2.uff'];

channel_data = uff.channel_data();
channel_data.read(filename,'/channel_data');

b_data_das = uff.beamformed_data();
b_data_cf = uff.beamformed_data();
b_data_pcf = uff.beamformed_data();
b_data_gcf = uff.beamformed_data();
b_data_mv = uff.beamformed_data();
b_data_ebmv = uff.beamformed_data();
b_data_dmas = uff.beamformed_data();
b_data_weights = uff.beamformed_data();

b_data_das.read(filename,'/b_data_das');
b_data_cf.read(filename,'/b_data_cf');
b_data_pcf.read(filename,'/b_data_pcf');
b_data_gcf.read(filename,'/b_data_gcf');
b_data_mv.read(filename,'/b_data_mv');
b_data_ebmv.read(filename,'/b_data_ebmv');
b_data_dmas.read(filename,'/b_data_dmas');
b_data_weights.read(filename,'/b_data_weights');

weights = b_data_weights.get_image('none');

%% Print authorship and citation details for the dataset
channel_data.print_authorship

%channel_data.name = {['New ',channel_data.name{1}]}

%channel_data.name = {['v4 ',channel_data.name{1}]}
%% DAS - Delay-And-Sum
% Create the individual images and the zoomed images on the boxes
mkdir([ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/'])


das_img = b_data_das.get_image('none').*weights;  % Compensation weighting
b_data_das.data = das_img(:); %Putting the weighted data back...
das_img = db(abs(das_img./max(das_img(:))));      % Normalize on max
f1 = figure(1);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,das_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/DAS'],'eps2c')
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/DAS'],'png')
f2 = figure(2);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,das_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
axis([-15 3 27 33]);
saveas(f2,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/DAS_zoomed'],'eps2c')


%%
image.all{1} = das_img;
image.all_signal{1} = double(b_data_das.get_image('none'));
image.all_signal{1} = (image.all_signal{1}./max(image.all_signal{1}(:)));
image.tags{1} = 'DAS';

z_start = 39;
z_stop = 48.5;
x_start = -20;
x_stop = 20;
[meanLines,~] = getMeanLate
ralLines(b_data_das,image,z_start,z_stop,-20,20);



figure(222);clf
subplot(211);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,b_data_das.get_image());colormap gray;caxis([-60 0]);
axis image;colorbar;ylim([39 48.5])
subplot(212);
plot(b_data_das.scan.x_axis*1000,meanLines.all{1},'LineWidth',2);hold all;
plot([-14 -14],[-80 0],'LineWidth',2);
plot([14 14],[-80 0],'LineWidth',2);
title('Average Response through full gradient');
%% GLT - Gray Level Transform
clear glt;
glt = postprocess.polynomial_gray_level_transform();
% glt.a = 4.3750e-04;
% glt.b = -0.0062;
% glt.c = 0.0500;
% glt.d = 0;
% glt.a = 0.500e-04;
% glt.b = -0.020;
% glt.c = 0.100;
% glt.d = 0;
% glt.a = 0.1e-04;
% glt.b = -0.025;
% glt.c = 0.100;
% glt.d = 0;

glt.a = 0.01e-04;
glt.b = -0.03;
glt.c = 0.100;
glt.d = 0;

% glt.a = 0.01e-04;
% glt.b = -0.011;
% glt.c = 0.500;
% glt.d = 0;

glt.plot_functions = 1;

glt.input = b_data_das;
glt.scan = b_data_das.scan;
b_data_glt = glt.go();


glt_img = b_data_glt.get_image('none');  % Compensation weighting
b_data_glt.data = glt_img(:);
glt_img = db(abs(glt_img./max(glt_img(:)))+eps);      % Normalize on max
f100 = figure(100);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,glt_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f100,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/GLT'],'eps2c')
saveas(f100,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/GLT'],'png')
f101 = figure(101);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,glt_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
axis([-15 3 27 33]);
saveas(f101,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/GLT_zoomed'],'eps2c')

%% Sigmoid curve

glt_s = postprocess.scurve_gray_level_transform();
% 
% glt_s.a = 0.08;
% glt_s.b = -35;
% glt_s.c = 0.008;

% glt_s.a = 0.08;
% glt_s.b = -35;
% glt_s.c = 0.008;

glt_s.a = 0.12;
glt_s.b = -40;
glt_s.c = 0.008;

% glt_s.a = 0.15;
% glt_s.b = -32.5;
% glt_s.c = 0.01;

glt_s.plot_functions = 1;
glt_s.input = b_data_das;
glt_s.scan = b_data_das.scan;
b_data_glt = glt_s.go();

%%

glt_img = b_data_glt.get_image('none');  % Compensation weighting
b_data_glt.data = glt_img(:);
glt_img = db(abs(glt_img./max(glt_img(:)))+eps);      % Normalize on max
f100 = figure(100);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,glt_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f100,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/GLT'],'eps2c')
saveas(f100,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/GLT'],'png')
f101 = figure(101);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,glt_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
axis([-15 3 27 33]);
saveas(f101,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/GLT_zoomed'],'eps2c')
%%
diff{8} = tools.dynamic_range_test(channel_data,b_data_glt,'GLT')

%% CF - Coherence Factor
cf_img = b_data_cf.get_image('none').*weights;
b_data_cf.data = cf_img(:);
cf_img = db(abs(cf_img./max(cf_img(:))));                 % Normalize on max
f2 = figure(2);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,cf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f2,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/CF'],'eps2c')
f22 = figure(22);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,cf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
axis([-15 3 27 33]);
saveas(f22,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/CF_zoomed'],'eps2c')

%% PCF - Phase Coherence Factor
pcf_img = b_data_pcf.get_image('none').*weights;
b_data_pcf.data = pcf_img(:);
b_data_pcf.data(isinf(b_data_pcf.data)) = eps;
b_data_pcf.data(b_data_pcf.data==0) = eps;
pcf_img(isinf(pcf_img)) = eps;
pcf_img(pcf_img==0) = eps;
pcf_img = db(abs(pcf_img./max(pcf_img(:)))+eps);
f3 = figure(3);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,pcf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f3,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/PCF'],'eps2c')
f33 = figure(33);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,pcf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
axis([-15 3 27 33]);
saveas(f33,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/PCF_zoomed'],'eps2c')

%% GCF - Generalized Coherence Factor
gcf_img = b_data_gcf.get_image('none').*weights;
b_data_gcf.data = gcf_img(:);
gcf_img = db(abs(gcf_img./max(gcf_img(:))));
f4 = figure(4);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,gcf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f4,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/GCF'],'eps2c')
f44 = figure(44);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,gcf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
axis([-15 3 27 33]);
saveas(f44,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/GCF_zoomed'],'eps2c')

%% MV - Capon's Minimum Variance
mv_img = b_data_mv.get_image('none').*weights;
b_data_mv.data = mv_img(:);
mv_img = db(abs(mv_img./max(mv_img(:))));
f6 = figure(6);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,mv_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f6,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/MV'],'eps2c')
f66 = figure(66);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,mv_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
axis([-15 3 27 33]);
saveas(f66,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/MV_zoomed'],'eps2c')

%% F-DMAS - Filtered-Delay-Multiply-and-Sum
dmas_img = b_data_dmas.get_image('none').*weights;
b_data_dmas.data = dmas_img(:);
dmas_img = db(abs(dmas_img./max(dmas_img(:))));
f7 = figure(7);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,dmas_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f7,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/DMAS'],'eps2c')
f77 = figure(77);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,dmas_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
axis([-15 3 27 33]);
saveas(f77,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/DMAS_zoomed'],'eps2c')

%% EBMV - Eigenspace Based Minimum Variance
ebmv_img = b_data_ebmv.get_image('none').*weights;
b_data_ebmv.data = ebmv_img(:);
ebmv_img = db(abs(ebmv_img./max(ebmv_img(:))));
f8 = figure(8);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,ebmv_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f8,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/EBMV'],'eps2c')
f88 = figure(88);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,ebmv_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
axis([-15 3 27 33]);
saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/EBMV_zoomed'],'eps2c')


%% Putting the dB images and the signal in a cell array
addpath([ustb_path,'/publications/DynamicRage/functions']);
image.all{1} = das_img;
image.all_signal{1} = double(b_data_das.get_image('none'));
image.all_signal{1} = (image.all_signal{1}./max(image.all_signal{1}(:)));
image.tags{1} = 'DAS';
image.all{2} = mv_img;
image.all_signal{2} = double(b_data_mv.get_image('none'));
image.all_signal{2} = (image.all_signal{2}./max(image.all_signal{2}(:)));
image.tags{2} = 'MV';
image.all{3} = ebmv_img;
image.all_signal{3} = double(b_data_ebmv.get_image('none'));
image.all_signal{3} = (image.all_signal{3}./max(image.all_signal{3}(:)));
image.tags{3} = 'EBMV';
image.all{4} = dmas_img;
image.all_signal{4} = double(b_data_dmas.get_image('none'));
image.all_signal{4} = (image.all_signal{4}./max(image.all_signal{4}(:)));
image.tags{4} = 'F-DMAS';
image.all{5} = cf_img;
image.all_signal{5} = double(b_data_cf.get_image('none'));
image.all_signal{5} = (image.all_signal{5}./max(image.all_signal{5}(:)));
image.tags{5} = 'CF';
image.all{6} = gcf_img;
image.all_signal{6} = double(b_data_gcf.get_image('none'));
image.all_signal{6} = (image.all_signal{6}./max(image.all_signal{6}(:)));
image.tags{6} = 'GCF';
image.all{7} = pcf_img;
image.all_signal{7} = double(b_data_pcf.get_image('none'));
image.all_signal{7}(isinf(image.all_signal{7})) = eps;
image.all_signal{7} = (image.all_signal{7}./max(image.all_signal{7}(:)));
image.tags{7} = 'PCF';
image.all{8} = glt_img;
image.all_signal{8} = double(b_data_glt.get_image('none'));
image.all_signal{8} = (image.all_signal{8}./max(image.all_signal{8}(:)));
image.tags{8} = 'GLT';

%% Calibrate images

figure;
subplot(211)
imagesc(image.all{7})
subplot(212)
imagesc(isinf(image.all{7}))

%%

% Using "linspecer" colors
colors=    [0.9047    0.1918    0.1988; ...
            0.2941    0.5447    0.7494; ...
            0.3718    0.7176    0.3612; ...
            1.0000    0.5482    0.1000; ...
            0.8650    0.8110    0.4330; ...
            0.6859    0.4035    0.2412; ...
            0.9718    0.5553    0.7741; ...
            0.6400    0.6400    0.6400];
        
        %%
        
[image_calibrated_on_theoretical] = calibrateImage_on_theoretical(b_data_das,image,40,48.5,-14,14,3,1.8,colors);

%%

[image_calibrated] = calibrateImage(b_data_das,image,40,48.5,-14,14,3,1.8,colors);


%%
x_axis_grad_start = 14;
mask_lateral=abs(b_data_das.scan.x_axis)<x_axis_grad_start*10^-3;
theory_lateral = (-1.8*(b_data_das.scan.x_axis(mask_lateral)+x_axis_grad_start*10^-3))*10^3;


[meanLines_full_calibrated,~] = getMeanLateralLines(b_data_das,image_calibrated,39,48.5,-x_axis_grad_start,x_axis_grad_start);

f10 = figure(77);clf;hold all;
subplot(211);
plot(theory_lateral,theory_lateral,'k--','DisplayName','Theory','linewidth',2); hold on; grid on;
for i = 1:length(image_calibrated.all)
    plot(theory_lateral,meanLines_full_calibrated.all{i}-max(meanLines_full_calibrated.all{i}),'linewidth',2,'DisplayName',image.tags{i},'Color',colors(i,:));
end
xlabel('Theoretical [dB]');
ylabel('Output [dB]');
set(gca, 'XDir','reverse');%title(image.tags{i});
xlim([-50 0])
ylim([-80 0])
legend show
%%
saveas(f10,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/calibrated/all_curves'],'eps2c')

%%
x_axis_grad_start = 14;
mask_lateral=abs(b_data_das.scan.x_axis)<x_axis_grad_start*10^-3;
theory_lateral = (-1.8*(b_data_das.scan.x_axis(mask_lateral)+x_axis_grad_start*10^-3))*10^3;


[meanLines_full_calibrated,~] = getMeanLateralLines(b_data_das,image_calibrated_on_theoretical,39,48.5,-x_axis_grad_start,x_axis_grad_start);

f10 = figure(77);clf;hold all;
subplot(211);
plot(theory_lateral,theory_lateral,'k--','DisplayName','Theory','linewidth',2); hold on; grid on;
for i = 1:length(image_calibrated_on_theoretical.all)
    plot(theory_lateral,meanLines_full_calibrated.all{i}-max(meanLines_full_calibrated.all{i}),'linewidth',2,'DisplayName',image.tags{i},'Color',colors(i,:));
end
xlabel('Theoretical [dB]');
ylabel('Output [dB]');
set(gca, 'XDir','reverse');%title(image.tags{i});
xlim([-50 0])
ylim([-80 0])
legend show
%%
saveas(f10,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/calibrated/all_curves_on_theoretical'],'eps2c')

%%
dynamic_range = -60;
%%
 figure(56566);
for i = 1:length(image.all)
    i
    subplot(3,8,i);
    imagesc(image.all{i});colormap gray;caxis([dynamic_range 0]);title(image.tags{i});
    subplot(3,8,8+i);
    imagesc(image_calibrated.all{i});colormap gray;caxis([dynamic_range 0]);title(image_calibrated.tags{i});
    %subplot(3,8,16+i);
    %imagesc(db(image_calibrated.signal_all{i}));colormap gray;caxis([dynamic_range 0]);title(image_calibrated.tags{i});
end
%%
dynamic_range = -60;
for i = 1:length(image.all)
    f10 = figure(10);clf;
    imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,image_calibrated.all{i});
    colormap gray;caxis([dynamic_range 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
    title(image_calibrated.tags{i});
    set(gca,'FontSize',14)
    saveas(f10,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/calibrated/',image.tags{i}],'eps2c')
    axis([-15 3 27 33]);
    saveas(f10,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/calibrated/',image.tags{i},'_zoomed'],'eps2c')
end

%%
f10 = figure(10);clf;
for i = 1:length(image.all)
    subplot(2,4,i);
    imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,image_calibrated.all{i});
    colormap gray;caxis([dynamic_range 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
    title(image_calibrated.tags{i});
    set(gca,'FontSize',14)
    axis([-15 3 27 33]);
end
%%
saveas(f10,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/calibrated/calibrated_boxes'],'eps2c')

%%
lines = getMeanLateralLines(b_data_das,image,16,18,-15,5)

figure(12123);clf;hold all;
for i = 1:length(image.all)

    plot(lines.all{i},'DisplayName',image.tags{i},'LineWidth',2);
end
legend


%%

%% Plot the lateral and axial gradient together with the experimental lateral gradient
image_experimental = load('Experimental.mat','image');
image_experimental = image_experimental.image;
b_data_das_exp = load('Experimental.mat','b_data_das'); 
b_data_das_exp = b_data_das_exp.b_data_das;

x_axis_grad_start = 14;

[meanLinesLateralExp,~] = getMeanLateralLines(b_data_das_exp,image_experimental,39,48.5,-x_axis_grad_start,x_axis_grad_start);

gradient = -1.8;

mask_lateral_exp=abs(b_data_das_exp.scan.x_axis)<x_axis_grad_start*10^-3;
%theory_lateral_exp=-40*(b_data_das_exp.scan.x_axis(mask_lateral_exp)+x_axis_grad_start*10^-3)/24.1e-3;
theory_lateral_exp = (gradient*(b_data_das_exp.scan.x_axis(mask_lateral_exp)+x_axis_grad_start*10^-3))*10^3;

[meanLinesLateral,x_axis] = getMeanLateralLines(b_data_das,image,39,48.5,-x_axis_grad_start,x_axis_grad_start);
mask_lateral=abs(b_data_das.scan.x_axis)<x_axis_grad_start*10^-3;
theory_lateral = (gradient*(b_data_das.scan.x_axis(mask_lateral)+x_axis_grad_start*10^-3))*10^3;

[meanLines_axial,z_axis] = getMeanAxialLines(b_data_das,image,10,38,15,18.5);
mask_axial= b_data_das.scan.z_axis<38e-3 & b_data_das.scan.z_axis>10e-3;
theory_axial = (gradient*(b_data_das.scan.z_axis(mask_axial)-10e-3))*10^3;

for i = 1:length(image.all)
    f88 = figure(8888+i);clf;hold all;
    plot(theory_lateral,meanLinesLateral.all{i}-max(meanLinesLateral.all{i}),'LineStyle','-.','Linewidth',2,'DisplayName','Sim. lateral');
    plot(theory_axial,meanLines_axial.all{i}-max(meanLines_axial.all{i}),'Linestyle','--','Linewidth',2,'DisplayName','Sim. axial');
    plot(theory_lateral_exp,meanLinesLateralExp.all{i}-max(meanLinesLateralExp.all{i}),'Linestyle',':','Linewidth',2,'DisplayName','Exp. lateral');
    plot(theory_lateral,theory_lateral,'LineStyle','-','Linewidth',2,'Color','k','DisplayName','Theory');
    set(gca, 'XDir','reverse');%title(image.tags{i});
    axis image
    legend('show','Location','sw');
    ylim([-60 0])
    xlim([-50 0])
    grid on
    xlabel('Theoretical [dB]');
    ylabel('Output [dB]');
    saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/gradient/',image.tags{i}],'eps2c')
    %saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/gradient/',image.tags{i}],'png')
end

%% Plot the lateral and axial gradient together with the experimental lateral gradient after calibration
% x_axis_grad_start = 15;
% 
% [meanLinesLateral_calib,x_axis] = getMeanLateralLines(b_data_das,image_calibrated,39,48.5,-x_axis_grad_start,x_axis_grad_start);
% mask_lateral=abs(b_data_das.scan.x_axis)<x_axis_grad_start*10^-3;
% theory_lateral = (-1.8*(b_data_das.scan.x_axis(mask_lateral)+x_axis_grad_start*10^-3))*10^3;
% 
% [meanLines_axial_calib,z_axis] = getMeanAxialLines(b_data_das,image_calibrated,10,40,15,18.5);
% mask_axial= b_data_das.scan.z_axis<40e-3 & b_data_das.scan.z_axis>10e-3;
% theory_axial = (-1.8*(b_data_das.scan.z_axis(mask_axial)-10e-3))*10^3;
% f88 = figure(9888);clf;hold all;
% for i = 1:length(image.all)
%     subplot(2,4,i);hold all;
%     plot(theory_lateral,meanLinesLateral.all{i}-mean(meanLinesLateral.all{i}(1:30)),'LineStyle','-.','Linewidth',2,'DisplayName','Uncalibrated');
%     plot(theory_lateral,meanLinesLateral_calib.all{i}-mean(meanLinesLateral_calib.all{i}(1:30)),'LineStyle','-.','Linewidth',2,'DisplayName','Calibrated');
%     %plot(theory_axial,meanLines_axial_calib.all{i}-max(meanLines_axial_calib.all{i}),'Linestyle','--','Linewidth',2,'DisplayName','Calibrated axial');
%     plot(theory_lateral,theory_lateral,'LineStyle','-','Linewidth',2,'Color','k','DisplayName','Theory');
%     set(gca, 'XDir','reverse');title(image.tags{i});
%     axis image
%     legend('show','Location','sw');
%     ylim([-60 0])
%     xlim([-50 0])
%     grid on
%     xlabel('Theoretical [dB]');
%     ylabel('Output [dB]');
% end
% saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/gradient/calibrated'],'eps2c')
% 
% f88 = figure(9889);clf;hold all;
% for i = 1:length(image.all)
%     subplot(2,4,i);hold all;
%     plot(theory_lateral,meanLinesLateral.all{i}-mean(meanLinesLateral.all{i}(1:30)),'LineStyle','-.','Linewidth',2,'DisplayName','Uncalibrated');
%     plot(theory_lateral,meanLinesLateral_calib.all{i}-mean(meanLinesLateral_calib.all{i}(1:30)),'LineStyle','-.','Linewidth',2,'DisplayName','Calibrated');
%     plot(theory_axial,meanLines_axial_calib.all{i}-max(meanLines_axial_calib.all{i}),'Linestyle','--','Linewidth',2,'DisplayName','Calibrated axial');
%     plot(theory_lateral,theory_lateral,'LineStyle','-','Linewidth',2,'Color','k','DisplayName','Theory');
%     set(gca, 'XDir','reverse');title(image.tags{i});
%     axis image
%     legend('show','Location','sw');
%     ylim([-60 0])
%     xlim([-50 0])
%     grid on
%     xlabel('Theoretical [dB]');
%     ylabel('Output [dB]');
% end
% saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/gradient/calibrated_with_axial'],'eps2c')
%% Plot the lateral line through the boxes
[meanLines,x_axis] = getMeanLateralLines(b_data_das,image,27.5,32.5,b_data_das.scan.x(1)*10^3,b_data_das.scan.x(end)*10^3);
theoretical = [-100 ones(1,255)*0 ones(1,255)*-10 -100 ones(1,256)*-100 ones(1,256)*-100 -100 ones(1,255)*0 ones(1,255)*-35 -100];
x_axis_theoretical = linspace(-13,2,1536);

f89 = figure(89);clf;hold all;
subplot(211);hold on;
plot(x_axis_theoretical,theoretical,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(5,:));
plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(6,:));
plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(7,:));
plot(x_axis,meanLines.all{8}-max(meanLines.all{8}),'Linewidth',2,'DisplayName',image.tags{8},'Color',colors(8,:));
ylim([-80 0]);
xlim([-14 8]);
legend('Location','best'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');

%%
saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/boxes'],'eps2c')

%% Plot the lateral line through the boxes after calibration
[meanLines,x_axis] = getMeanLateralLines(b_data_das,image_calibrated,27.5,32.5,b_data_das.scan.x(1)*10^3,b_data_das.scan.x(end)*10^3);
theoretical = [-100 ones(1,255)*0 ones(1,255)*-10 -100 ones(1,256)*-100 ones(1,256)*-100 -100 ones(1,255)*0 ones(1,255)*-35 -100];
x_axis_theoretical = linspace(-13,2,1536);

f89 = figure(90);clf;hold all;
subplot(211);hold on;
plot(x_axis_theoretical,theoretical,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(5,:));
plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(6,:));
plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(7,:));
plot(x_axis,meanLines.all{8}-max(meanLines.all{8}),'Linewidth',2,'DisplayName',image.tags{8},'Color',colors(8,:));
ylim([-80 0]);
xlim([-14 8]);
legend('Location','best'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');

%%
saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/calibrated/boxes_calibrated_lines'],'eps2c')
saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/calibrated/boxes_calibrated_lines'],'png')
%% Plot the lateral line through the boxes after calibration
[meanLines,x_axis] = getMeanLateralLines(b_data_das,image_calibrated_on_theoretical,27.5,32.5,b_data_das.scan.x(1)*10^3,b_data_das.scan.x(end)*10^3);
theoretical = [-100 ones(1,255)*0 ones(1,255)*-10 -100 ones(1,256)*-100 ones(1,256)*-100 -100 ones(1,255)*0 ones(1,255)*-35 -100];
x_axis_theoretical = linspace(-13,2,1536);

f89 = figure(90);clf;hold all;
subplot(211);hold on;
plot(x_axis_theoretical,theoretical,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(5,:));
plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(6,:));
plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(7,:));
plot(x_axis,meanLines.all{8}-max(meanLines.all{8}),'Linewidth',2,'DisplayName',image.tags{8},'Color',colors(8,:));
ylim([-80 0]);
xlim([-14 8]);
legend('Location','best'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');

%%
saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/calibrated/boxes_calibrated_lines_on_theoretical'],'eps2c')
saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/calibrated/boxes_calibrated_lines_on_theoretical'],'png')

%% Measure contrast and plot together with the experimental
[CR_signal_sim, CR_signal_dagger, CR_image_sim, CNR_signal_sim, CNR_image_sim, CE_alfonso,CNRE_alfonso] = measureContrast(b_data_das,image,-5.5,17.5,3.5,5,8,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/DAS_cyst_indicated']);

load('Experimental.mat','CR_signal_exp','CNR_signal_exp','CR_image_exp','CNR_image_exp');

%%
f9 = figure;
subplot(211);
bar([10*log10(CR_signal_sim)' 10*log10(CR_signal_exp)'])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CR [dB]');
x_pos = [0.8 1.8 2.8 3.8 4.8 5.8 6.8 7.8];
text(x_pos,double(10*log10(CR_signal_sim(:))),num2str(round(10*log10(CR_signal_sim(:))),'%d'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
x_pos = [1.2 2.2 3.2 4.2 5.2 6.2 7.2 8.2];
text(x_pos,double(10*log10(CR_signal_exp(:))),num2str(round(10*log10(CR_signal_exp(:))),'%d'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
set(gca,'Ydir','reverse')
legend('Sim.','Exp.')
grid on

%%
saveas(f9,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/CR'],'eps2c')
%%
f90 = figure;
subplot(212);
bar([CR_image_sim' CR_image_exp'])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CR_{LC} [dB]');
x_pos = [0.8 1.8 2.8 3.8 4.8 5.8 6.8 7.8];
text(x_pos,double(CR_image_sim(:)),num2str(round(CR_image_sim(:)),'%d'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
x_pos = [1.2 2.2 3.2 4.2 5.2 6.2 7.2 8.2];
text(x_pos,double(CR_image_exp(:)),num2str(round(CR_image_exp(:)),'%d'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
ylim([0 125])
legend('Sim.','Exp.')
grid on
%%
saveas(f90,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/CR_LC'],'eps2c')
%%
f77 = figure(77);clf;
subplot(211);
bar([CNR_signal_sim' CNR_signal_exp'])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CNR');
x_pos = [0.75 1.75 2.75 3.75 4.75 5.75 6.75 7.75];
text(x_pos,double(CNR_signal_sim(:)),num2str((CNR_signal_sim(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
x_pos = [1.25 2.25 3.25 4.25 5.25 6.25 7.25 8.25];
text(x_pos,double(CNR_signal_exp(:)),num2str((CNR_signal_exp(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
ylim([0 1.45])
legend('Sim.','Exp.')
grid on

%%
saveas(f77,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/CNR'],'eps2c')
%%
f78 = figure(78);clf;
subplot(212);
bar([CNR_image_sim' CNR_image_exp'])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CNR_{LC}');
x_pos = [0.75 1.75 2.75 3.75 4.75 5.75 6.75 7.75];
text(x_pos,double(CNR_image_sim(:)),num2str((CNR_image_sim(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
x_pos = [1.25 2.25 3.25 4.25 5.25 6.25 7.25 8.25];
text(x_pos,double(CNR_image_exp(:)),num2str((CNR_image_exp(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
legend('Sim.','Exp.')
grid on
ylim([0 9])
%%

saveas(f78,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/CNR_LC'],'eps2c')

%% Find the corrlation between dynamic range stretching and measured contrast as CR
figure(77);clf;
for i = 1:length(image.all)    
    simulated_lateral = meanLinesLateral.all{i}-max(meanLinesLateral.all{i});
    experimental_lateral = meanLinesLateralExp.all{i}-max(meanLinesLateralExp.all{i});
    
    % Remove any NAN's (we have one for PCF)
    experimental_lateral = experimental_lateral(find(~isinf(experimental_lateral)));
    theory_lateral_exp_temp = theory_lateral_exp(find(~isinf(experimental_lateral))); 
    
    % Do a polynomial fit to get an estimate of the gradient
    [sim_poly_fit,S_sim,mu_sim] = polyfit(theory_lateral, simulated_lateral',3);
    [exp_poly_fit,S_exp,mu_exp] = polyfit(theory_lateral_exp_temp, experimental_lateral',3);
   
    
    x = -40; %Calculate devation for -40 dB
    y_sim(i) = polyval(sim_poly_fit,x,S_sim,mu_sim)-x;
    y_exp(i) = polyval(exp_poly_fit,x,S_exp,mu_exp)-x;
        
    % Calculate difference in CR value compared to DAS
    diff_CR_vs_das_sim(i) = CR_image_sim(i) - CR_image_sim(1);
    diff_CR_vs_das_exp(i) = CR_image_exp(i) - CR_image_exp(1);
    
    x_vector = linspace(0,-40,length(simulated_lateral));
    %Plot the fitted polynomial to see the match
    figure(77);hold all;
    subplot(211);hold all;
    plot(simulated_lateral,'DisplayName',[image.tags{i}, ' Sim. gradient'],'Color',colors(i,:));
    plot(polyval(sim_poly_fit,theory_lateral,S_sim,mu_sim),'DisplayName',[image.tags{i},' fitted polynomial'],'Color',colors(i,:));
    title('Simulated');
    subplot(212);hold all;
    plot(experimental_lateral,'DisplayName',[image.tags{i},' Exp. gradient'],'Color',colors(i,:));
    plot(polyval(exp_poly_fit,theory_lateral_exp_temp,S_exp,mu_exp),'DisplayName',[image.tags{i}, ' fitted polynomial'],'Color',colors(i,:));
    title('Experimental');
end
set(gcf,'Position',[179 113 782 549]);

%% Create correlation plot for the simulated data
mdl_sim = fitlm(diff_CR_vs_das_sim(1:8),y_sim(1:8));
mdl_exp = fitlm(diff_CR_vs_das_exp(1:8),y_exp(1:8));

rsquared_ord_sim = mdl_sim.Rsquared.Ordinary
rsquared_adj_sim = mdl_sim.Rsquared.Adjusted

f111 = figure(111);clf;hold all;
linespes = '+ox*sd<>';
for i = 1:length(image.all)
    plot(diff_CR_vs_das_sim(i),y_sim(i),linespes(i),'Color',colors(i,:),'MarkerSize',10);
end
mdl_sim.plot
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',3)
set(hline(3),'Color',[0 0 0]);
set(hline(3),'DisplayName','Fitted line');
title('');
xlabel('CR diff. from DAS');xlim([-5 110]);
ylabel('Deviation at -40 dB');ylim([-85 5]);
text(40,0,sprintf('R-Squared: %.2f',rsquared_ord_sim),'FontSize',18);
text(40,-5,sprintf('R-Squared adj: %.2f',rsquared_adj_sim),'FontSize',18);
set(gca,'FontSize',14);
delete(hline(4));
legend();
hLegend = findobj(gcf, 'Type', 'Legend');
hLegend.String;
legend show;
legend([image.tags hLegend.String],'Location','sw');

saveas(f111,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/correlation_sim'],'eps2c');

%% Create correlation plot for the experimental data
rsquared_ord_exp = mdl_exp.Rsquared.Ordinary
rsquared_adj_exp = mdl_exp.Rsquared.Adjusted

f112 = figure(112);clf;hold all;
linespes = '+ox*sd<>';
for i = 1:length(image.all)
    plot(diff_CR_vs_das_exp(i),y_exp(i),linespes(i),'Color',colors(i,:),'MarkerSize',10);
end
mdl_exp.plot
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',3);
set(hline(3),'Color',[0 0 0]);
set(hline(3),'DisplayName','Fitted line');
title('')
xlabel('CR diff. from DAS');xlim([-5 75]);
ylabel('Deviation at -40 dB');ylim([-50 5]);
text(40,0,sprintf('R-Squared: %.2f',rsquared_ord_exp),'FontSize',18);
text(40,-3,sprintf('R-Squared adj: %.2f',rsquared_adj_exp),'FontSize',18);
set(gca,'FontSize',14);
delete(hline(4));
legend();
hLegend = findobj(gcf, 'Type', 'Legend');
legend show
legend([image.tags hLegend.String],'Location','sw');

saveas(f112,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/correlation_exp'],'eps2c');

%%




diff{1} = tools.dynamic_range_test(channel_data,b_data_das,'DAS')
set(gcf,'Position',[7 222 1271 347])
saveas(gcf,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/dynamic_range_test/','DAS'],'eps2c')
diff{2} = tools.dynamic_range_test(channel_data,b_data_mv,'MV')
set(gcf,'Position',[7 222 1271 347])
saveas(gcf,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/dynamic_range_test/','MV'],'eps2c')
diff{3} = tools.dynamic_range_test(channel_data,b_data_ebmv,'EBMV')
set(gcf,'Position',[7 222 1271 347])
saveas(gcf,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/dynamic_range_test/','EBMV'],'eps2c')
diff{4} = tools.dynamic_range_test(channel_data,b_data_dmas,'F-DMAS')
set(gcf,'Position',[7 222 1271 347])
saveas(gcf,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/dynamic_range_test/','F-DMAS'],'eps2c')
diff{5} = tools.dynamic_range_test(channel_data,b_data_cf,'CF')
set(gcf,'Position',[7 222 1271 347])
saveas(gcf,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/dynamic_range_test/','CF'],'eps2c')
diff{6} = tools.dynamic_range_test(channel_data,b_data_gcf,'GCF')
set(gcf,'Position',[7 222 1271 347])
saveas(gcf,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/dynamic_range_test/','GCF'],'eps2c')
diff{7} = tools.dynamic_range_test(channel_data,b_data_pcf,'PCF')
set(gcf,'Position',[7 222 1271 347])
saveas(gcf,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/dynamic_range_test/','PCF'],'eps2c')
diff{8} = tools.dynamic_range_test(channel_data,b_data_glt,'GLT')
set(gcf,'Position',[7 222 1271 347])
saveas(gcf,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/dynamic_range_test/','GLT'],'eps2c')

%%
for i = 1:length(image.tags)
    diff_avrg(i) = (mean(diff{i}));
    diff_lateral(i) = diff{i}(1);
    diff_axial(i) = diff{i}(2);
    diff_CR_vs_das_sim(i) = CR_image_sim(i) - CR_image_sim(1);
    diff_CR_signal_vs_das_sim(i) = CR_signal_sim(i) - CR_signal_sim(1);
end

%%
cr_percentage = (abs((CR_image_sim - CR_image_sim(1))/CR_image_sim(1)))*100;
cr_percentage = CR_image_sim - CR_image_sim(1);
%cr_percentage = (abs((10*log10(CR_signal_sim) - 10*log10(CR_signal_sim(1)))/10*log10(CR_signal_sim(1))))*100;

CR_selected = cr_percentage;%abs(diff_CR_vs_das_sim);

lmd_sim_slope = fitlm(CR_selected(1:8),diff_avrg(1:8))

lmd_sim_slope.Rsquared.Ordinary
lmd_sim_slope.Rsquared.Adjusted

f113 = figure(114);clf;hold all;
linespes = '+ox*sd<>';
for i = 1:length(image.all)
    plot(CR_selected(i),diff_avrg(i),linespes(i),'Color',colors(i,:),'MarkerSize',10);
end
lmd_sim_slope.plot
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',3);
set(hline(3),'Color',[0 0 0]);
set(hline(3),'DisplayName','Fitted line');
title('')
xlabel('CR improve cmp. DAS [%]');%xlim([-5 75]);
ylabel('Slope diff [in %]');%ylim([-50 5]);
text(80,10,sprintf('R-Squared: %.2f',lmd_sim_slope.Rsquared.Ordinary),'FontSize',18);
text(80,0,sprintf('R-Squared adj: %.2f',lmd_sim_slope.Rsquared.Adjusted),'FontSize',18);
set(gca,'FontSize',14);
delete(hline(4));
legend();
hLegend = findobj(gcf, 'Type', 'Legend');
legend show
legend([{image.tags{1:8}} hLegend.String],'Location','nw');
figure;
lmd_sim_slope.plot();

saveas(f113,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/correlation_simulation'],'eps2c')
%%
lmd_sim_slope = fitlm(CR_selected(1:7),diff_avrg(1:7))

lmd_sim_slope.Rsquared.Ordinary
lmd_sim_slope.Rsquared.Adjusted

f113 = figure(114);clf;hold all;
linespes = '+ox*sd<>';
for i = 1:length(image.all)-1
    plot(CR_selected(i),diff_avrg(i),linespes(i),'Color',colors(i,:),'MarkerSize',10);
end
lmd_sim_slope.plot
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',3);
set(hline(3),'Color',[0 0 0]);
set(hline(3),'DisplayName','Fitted line');
title('')
xlabel('CR_{LC} improve compared to DAS [dB]','Interpreter', 'tex');
ylabel('DRT value');
text(40,1,sprintf('R-Squared: %.2f',lmd_sim_slope.Rsquared.Ordinary),'FontSize',18);
text(40,0.9,sprintf('R-Squared adj: %.2f',lmd_sim_slope.Rsquared.Adjusted),'FontSize',18);
set(gca,'FontSize',14);
delete(hline(4));
legend();
hLegend = findobj(gcf, 'Type', 'Legend');
legend show
legend([{image.tags{1:7}} hLegend.String],'Location','nw');


saveas(f113,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/correlation_simulation_without_glt'],'eps2c')

%%
load('Experimental.mat','diff_avrg_exp');
f10 = figure(10);clf;
subplot(211);
bar([diff_avrg' diff_avrg_exp'])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('DRT value');
legend('Location','nw','sim.','exp.')
x_pos = [0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75 7.25 7.75 8.25];
drt_diff_all = [];
for i = 1:length(diff_avrg_exp) 
    drt_diff_all = [drt_diff_all; diff_avrg(i); diff_avrg_exp(i)];
end
text(x_pos,double(drt_diff_all),num2str(drt_diff_all(:),'%0.1f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
grid on;
ylim([0 2.7])

saveas(f10,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/DRT'],'eps2c')

%%
diff_CR_vs_das_sim_no_ebmv = [diff_CR_vs_das_sim(1:2) diff_CR_vs_das_sim(3:7)];
diff_lateral_no_ebmv = [diff_avrg(1:2) diff_avrg(3:7)];

lmd_sim_slope = fitlm(diff_CR_vs_das_sim_no_ebmv,diff_lateral_no_ebmv)

lmd_sim_slope.Rsquared.Ordinary
lmd_sim_slope.Rsquared.Adjusted

%%
%[handle_1] = plotLateralLine(b_data_tx,image,15,-9,-6,colors)
[handle_2] = plotLateralLine(b_data_das,image,35,-7.5,-3.5,colors)

fwhm_1 = measureFWHM(b_data_das,image,35,-7.5,-3.5)
fwhm_2 = measureFWHM(b_data_das,image,35,-7.5,-3.5)


f9 = figure(9);clf
ax = subplot(2,1,1);
% copyobj(allchild(handle_1),ax)
% ylim([-80 0]);xlim([-9 -6]); 
% legend('show','Location','best');
% ylabel('Amplitude [dB]');xlabel('x [mm]');
% title('PSF at z = 15 mm');
ax = subplot(2,1,1);
copyobj(allchild(handle_2),ax)
ylim([-80 0]);xlim([-7.5,-3.5]); 
legend('show','Location','best');
ylabel('Amplitude [dB]');xlabel('x [mm]');
grid on
%title('PSF at z = 45 mm');
%%
saveas(f9,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/PSF_at_z35'],'eps2c')

%%
f10 = figure(10);clf;
subplot(211);
bar([fwhm_1' fwhm_2'])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('FWHM [mm]');
legend('Location','best','PSF at z = 15 mm','PSF at z = 45 mm')
x_pos = [0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75 7.25 7.75 8.25];
fwhm_all = [];
for i = 1:length(fwhm_1) 
    fwhm_all = [fwhm_all; fwhm_1(i); fwhm_2(i)];
end
text(x_pos,double(fwhm_all(:)),num2str(fwhm_all(:),'%0.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
grid on;
%%
saveas(f10,[ustb_path,filesep,'publications/DynamicRage/figures/simulation_v2/FWHM'],'eps2c')
%%
selection = [1 4 5 6 7];
%selection = 1:8;
fwhm_2_alt = fwhm_2(selection);
fwhm_2_diff = fwhm_2_alt-fwhm_2_alt(1);
diff_avrg_alt = diff_avrg(selection);



lmd_sim_fwhm = fitlm(fwhm_2_diff,diff_avrg_alt)

lmd_sim_fwhm.Rsquared.Ordinary
lmd_sim_fwhm.Rsquared.Adjusted

f114 = figure(114);clf;hold all;
linespes = '+ox*sd<>';
for i = 1:length(fwhm_2_diff)
    plot(fwhm_2_diff(i),diff_avrg(i),linespes(i),'Color',colors(i,:),'MarkerSize',10);
end
lmd_sim_fwhm.plot
hline = findobj(gcf, 'type', 'line');
set(hline,'LineWidth',3);
set(hline(3),'Color',[0 0 0]);
set(hline(3),'DisplayName','Fitted line');
title('')
xlabel('CR diff. from DAS');%xlim([-5 75]);
ylabel('Slope diff [in %]');%ylim([-50 5]);
text(0.35,400,sprintf('R-Squared: %.2f',lmd_sim_fwhm.Rsquared.Ordinary),'FontSize',18);
text(0.35,350,sprintf('R-Squared adj: %.2f',lmd_sim_fwhm.Rsquared.Adjusted),'FontSize',18);
set(gca,'FontSize',14);
delete(hline(4));
legend();
hLegend = findobj(gcf, 'Type', 'Legend');
legend show
legend([image.tags{selection} hLegend.String],'Location','sw');
%lmd_sim_fwhm.plot();