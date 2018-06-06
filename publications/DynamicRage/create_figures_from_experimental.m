%% Create the figures from the experimental dynamic range phantom
% For the publication ?Rindal, O. M. H., Austeng, A., Fatemi, A., 
% & Rodriguez-Molares, A. (2018). The effect of dynamic range transformations
% in the estimation of contrast. Submitted to IEEE Transactions on Ultrasonics,
% Ferroelectrics, and Frequency Control.

clear all;
close all;
%% Load the data
filename = [data_path,filesep,'experimental_dynamic_range_phantom.uff'];

channel_data = uff.channel_data();
channel_data.read(filename,'/channel_data')

b_data_tx = uff.beamformed_data();
b_data_das = uff.beamformed_data();
b_data_cf = uff.beamformed_data();
b_data_pcf = uff.beamformed_data();
b_data_gcf = uff.beamformed_data();
b_data_mv = uff.beamformed_data();
b_data_ebmv = uff.beamformed_data();
b_data_dmas = uff.beamformed_data();
b_data_glt = uff.beamformed_data();

b_data_tx.read(filename,'/b_data_tx');
b_data_das.read(filename,'/b_data_das');
b_data_cf.read(filename,'/b_data_cf');
b_data_pcf.read(filename,'/b_data_pcf');
b_data_gcf.read(filename,'/b_data_gcf');
b_data_mv.read(filename,'/b_data_mv');
b_data_ebmv.read(filename,'/b_data_ebmv');
b_data_dmas.read(filename,'/b_data_dmas');
b_data_glt.read(filename,'/b_data_glt');

%% Print authorship and citation details for the dataset
print_authorship

%% Display and save the images from all beamformers under study
mkdir([ustb_path,filesep,'publications/DynamicRage/figures/experimental/'])

%% DAS
das_img = b_data_das.get_image('none'); 
das_img = db(abs(das_img./max(das_img(:)))); % Normalize on max
f1 = figure(1);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,das_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/DAS'],'eps2c')
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/DAS'],'png')

%% CF
cf_img = b_data_cf.get_image('none');
cf_img = db(abs(cf_img./max(cf_img(:)))); % Normalize on max
f2 = figure(2);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,cf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f2,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/CF'],'eps2c')

%% PCF
pcf_img = b_data_pcf.get_image('none');
pcf_img = db(abs(pcf_img./max(pcf_img(:)))); % Normalize on max
f3 = figure(3);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,pcf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f3,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/PCF'],'eps2c')

%% GCF
gcf_img = b_data_gcf.get_image('none');
gcf_img = db(abs(gcf_img./max(gcf_img(:)))); % Normalize on max
f4 = figure(4);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,gcf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f4,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/GCF'],'eps2c')

%% MV
mv_img = b_data_mv.get_image('none');
mv_img = db(abs(mv_img./max(mv_img(:)))); % Normalize on max
f5 = figure(6);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,mv_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f5,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/MV'],'eps2c')

%% F-DMAS
dmas_img = b_data_dmas.get_image('none');
dmas_img = db(abs(dmas_img./max(dmas_img(:)))); % Normalize on max
f7 = figure(7);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,dmas_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f7,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/DMAS'],'eps2c')

%% EBMV
ebmv_img = b_data_ebmv.get_image('none');
ebmv_img = db(abs(ebmv_img./max(ebmv_img(:)))); % Normalize on max
f6 = figure(8);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,ebmv_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f6,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/EBMV'],'eps2c')

%% GLT
glt = postprocess.gray_level_transform();
glt.a = 4.3750e-04;
glt.b = -0.0062;
glt.c = 0.0500;
glt.d = 0;
glt.plot_functions = 1;

glt.input = b_data_das;
glt.scan = b_data_tx.scan;
b_data_glt = glt.go();

glt_img = b_data_glt.get_image('none');
glt_img = db(abs(glt_img./max(glt_img(:))));
f9 = figure(9);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,glt_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f9,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/GLT'],'eps2c')

%% Put the images in a cell array.
% Put both the dB image and the signal before envelope detection
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
image.all_signal{7} = (image.all_signal{7}./max(image.all_signal{7}(:)));
image.tags{7} = 'PCF';
image.all{8} = glt_img;
image.all_signal{8} = double(b_data_glt.get_image('none'));
image.all_signal{8} = (image.all_signal{8}./max(image.all_signal{8}(:)));
image.tags{8} = 'GLT';

%% Add the path to the functions used in the evaluation
addpath functions

% Add the colors to be used using "linspecer" colors
colors=    [0.9047    0.1918    0.1988; ...
            0.2941    0.5447    0.7494; ...
            0.3718    0.7176    0.3612; ...
            1.0000    0.5482    0.1000; ...
            0.8650    0.8110    0.4330; ...
            0.6859    0.4035    0.2412; ...
            0.9718    0.5553    0.7741; ...
            0.6400    0.6400    0.6400];



%% Plot lateral graident in db vs. db plot. 
% NB! The final figure is plotted together with the simulated data in the
%  create_figures_from_simulation.m script. We are just saving these values
%  for now.

[meanLinesLateral,x_axis] = getMeanLateralLines(b_data_das,image,39,48.5,-12.5,12.5);
mask_lateral=abs(b_data_tx.scan.x_axis)<12.5e-3;
theory_lateral=-40*(b_data_tx.scan.x_axis(mask_lateral)+12.5e-3)/25e-3;

mkdir([ustb_path,filesep,'publications/DynamicRage/figures/experimental/gradient/'])

for i = 1:length(image.all)
    f88 = figure(8888+i);clf;hold all;
    plot(theory_lateral,meanLinesLateral.all{i}-max(meanLinesLateral.all{i}),...
            'LineStyle','-.','Linewidth',2,'Color','b','DisplayName','Lateral');
    plot(theory_lateral,theory_lateral,'LineStyle','-','Linewidth',2,'Color',...
                                                    'k','DisplayName','Theory');
    set(gca, 'XDir','reverse');title(image.tags{i});
    axis image
    legend show
    ylim([-100 0])
    xlim([-50 0])
    grid on
    xlabel('Theoretical [dB]');
    ylabel('Output [dB]');
    saveas(f88,[ustb_path,filesep,...
    'publications/DynamicRage/figures/experimental/gradient/',image.tags{i}],'eps2c')
end

%% Measure contrast
% Once again we are just saving the values to make the final plot combined
% with the simulated data in the create_figures_from_simulation.m
[CR_signal_exp, CR_signal_dagger, CR_image_exp, CNR_signal_exp, CNR_image_exp] ...
    = measureContrast(b_data_tx,image,-5.5,17.5,3.5,5,8,[ustb_path,filesep,...
    'publications/DynamicRage/figures/experimental/DAS_cyst_indicated']);

f9 = figure;
subplot(211);
bar([10*log10(CR_signal_exp)])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CR_{signal} [dB]');
x_pos = [1 2 3 4 5 6 7 8];
text(x_pos,double(10*log10(CR_signal_exp(:))),num2str(round(10*log10(CR_signal_exp(:))),'%d'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
set(gca,'Ydir','reverse')

subplot(212);
bar([CR_image_exp])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CR_{image} [dB]');
x_pos = [1 2 3 4 5 6 7 8];
text(x_pos,double(CR_image_exp(:)),num2str(round(CR_image_exp(:)),'%d'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)

saveas(f9,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/CR'],'eps2c')

f77 = figure(77);
subplot(211);
bar([CNR_signal_exp])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CNR_{signal}');
x_pos = [1 2 3 4 5 6 7 8];
text(x_pos,double(CNR_signal_exp(:)),num2str((CNR_signal_exp(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
ylim([0 1.1])

subplot(212);
bar([CNR_image_exp])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CNR_{image}');
x_pos = [1 2 3 4 5 6 7 8];
text(x_pos,double(CNR_image_exp(:)),num2str((CNR_image_exp(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)

saveas(f77,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/CNR'],'eps2c')

save('Experimental.mat','b_data_das','image','CR_signal_exp', ...
                        'CR_image_exp', 'CNR_signal_exp', 'CNR_image_exp');