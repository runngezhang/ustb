clear all;
close all;

filename = [data_path,filesep,'experimental_gradient_phantom.uff'];
filename = '/Users/omrindal/Development/USTB_dynamic_range/data/experimental_gradient_phantom.uff'

ustb_path = '/Users/omrindal/Development/USTB_dynamic_range'

b_data_tx = uff.beamformed_data();
b_data_das = uff.beamformed_data();
b_data_cf = uff.beamformed_data();
b_data_pcf = uff.beamformed_data();
b_data_gcf = uff.beamformed_data();
b_data_mv = uff.beamformed_data();
b_data_ebmv = uff.beamformed_data();
b_data_dmas = uff.beamformed_data();
b_data_weights = uff.beamformed_data();

b_data_tx.read(filename,'/b_data_tx');
b_data_das.read(filename,'/b_data_das');
b_data_cf.read(filename,'/b_data_cf');
b_data_pcf.read(filename,'/b_data_pcf');
b_data_gcf.read(filename,'/b_data_gcf');
b_data_mv.read(filename,'/b_data_mv');
b_data_ebmv.read(filename,'/b_data_ebmv');
b_data_dmas.read(filename,'/b_data_dmas');
b_data_weights.read(filename,'/b_data_weights');

weights = b_data_weights.get_image('none');
%%

das_img = b_data_das.get_image('none').*weights;  % Compensation weighting
das_img = db(abs(das_img./max(das_img(:))));                 % Normalize on max
f1 = figure(1);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,das_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/DAS'],'eps2c')
%saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/DAS'],'png')

cf_img = b_data_cf.get_image('none').*weights;
cf_img = db(abs(cf_img./max(cf_img(:))));                 % Normalize on max
f2 = figure(2);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,cf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f2,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/CF'],'eps2c')

pcf_img = b_data_pcf.get_image('none').*weights;
pcf_img = db(abs(pcf_img./max(pcf_img(:))));
f3 = figure(3);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,pcf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f3,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/PCF'],'eps2c')

gcf_img = b_data_gcf.get_image('none').*weights;
gcf_img = db(abs(gcf_img./max(gcf_img(:))));
f4 = figure(4);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,gcf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f4,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/GCF'],'eps2c')

mv_img = b_data_mv.get_image('none').*weights;
mv_img = db(abs(mv_img./max(mv_img(:))));
f5 = figure(6);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,mv_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f5,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/MV'],'eps2c')

dmas_img = b_data_dmas.get_image('none').*weights;
dmas_img = db(abs(dmas_img./max(dmas_img(:))));
f7 = figure(7);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,dmas_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f7,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/DMAS'],'eps2c')

ebmv_img = b_data_ebmv.get_image('none').*weights;
ebmv_img = db(abs(ebmv_img./max(ebmv_img(:))));
f6 = figure(8);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,ebmv_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f6,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/EBMV'],'eps2c')


%%
image.all{1} = das_img;
image.tags{1} = 'DAS';
image.all{5} = mv_img;
image.tags{5} = 'MV';
image.all{6} = ebmv_img;
image.tags{6} = 'EBMV';
image.all{2} = cf_img;
image.tags{2} = 'CF';
image.all{7} = pcf_img;
image.tags{7} = 'PCF';
image.all{3} = gcf_img;
image.tags{3} = 'GCF';
image.all{4} = dmas_img;
image.tags{4} = 'DMAS';

%%
addpath functions
addpath cbrewer

load colorbrewer.mat

colors = colorbrewer.qual.Dark2{7}./255;
[meanLines,x_axis] = getMeanLateralLines(b_data_das,image,39,48.5,-12.5,12.5);

mask=abs(b_data_tx.scan.x_axis)<12.5e-3;
theory=-40*(b_data_tx.scan.x_axis(mask)+12.5e-3)/25e-3;

f89 = figure(89);clf
subplot(211);hold all;
plot(b_data_tx.scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));

ylim([-80 0]);
xlim([-12.5 12.5])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
%title('Lateral gradient');
set(gca,'Fontsize',14)

subplot(212);hold all;
plot(b_data_tx.scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(5,:));
plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(6,:));
plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(7,:));
ylim([-80 0]);
xlim([-12.5 12.5])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
%title('Lateral gradient');

set(gca,'Fontsize',14)
saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/gradient'],'eps2c')
%saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/gradient'],'png')

%%
[meanLines,x_axis] = getMeanLateralLines(b_data_das,image,39,48.5,-20,20);
mask=abs(b_data_tx.scan.x_axis)<12.5e-3;
theory=-7.5-40*(b_data_tx.scan.x_axis(mask)+12.5e-3)/25e-3;

f89 = figure(90);clf;hold all;
plot(b_data_tx.scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0 ]);
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(5,:));
plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(6,:));
plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(7,:));
ylim([-80 0]);
xlim([-20 20])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
%title('Lateral gradient');

set(gca,'Fontsize',20)
saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/gradient_2'],'eps2c')

%%
%[handle_1] = plotLateralLine(b_data_tx,image,10,10,14,colors)
%[handle_2] = plotLateralLine(b_data_tx,image,20,10,14,colors)
[handle_3] = plotLateralLine(b_data_tx,image,30,10,14,colors)
fwhm_1 = measureFWHM(b_data_tx,image,10,10,14)
fwhm_2 = measureFWHM(b_data_tx,image,20,10,14)
fwhm_3 = measureFWHM(b_data_tx,image,30,10,14)

%%
f9 = figure(9);clf
% ax = subplot(3,1,1);
% copyobj(allchild(handle_1),ax)
% ylim([-80 0]);xlim([10 14]); 
% legend('show','Location','best');
% ylabel('Amplitude [dB]');xlabel('x [mm]');
% title('PSF at z = 10 mm');
% ax = subplot(3,1,2);
% copyobj(allchild(handle_2),ax)
% ylim([-80 0]);xlim([10 14]); 
% legend('show','Location','best');
% ylabel('Amplitude [dB]');xlabel('x [mm]');
% title('PSF at z = 20 mm');
ax = subplot(2,1,1);
copyobj(allchild(handle_3),ax)
ylim([-80 0]);xlim([10 14]); 
legend('show','Location','best');
ylabel('Amplitude [dB]');xlabel('x [mm]');
%title('PSF at z = 30 mm');
saveas(f9,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/PSF'],'eps2c')

%%
f10 = figure(100);clf;
subplot(211);
bar([fwhm_1' fwhm_2' fwhm_3'])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
x_pos = [0.65 1 1.35 1.65 2 2.35  2.65 3 3.35 3.65 4 4.35 4.65 5 5.35 5.65 6 6.35 6.65 7 7.35];
fwhm_all = [];
for i = 1:length(fwhm_1) 
    fwhm_all = [fwhm_all; fwhm_1(i); fwhm_2(i); fwhm_3(i)];
end
text(x_pos,double(fwhm_all(:)),num2str(fwhm_all(:),'%0.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',10)
grid on;

ylabel('FWHM [mm]');
legend('PSF at z = 10 mm','PSF at z = 20 mm','PSF at z = 30 mm')
saveas(f10,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/FWHM'],'eps2c')
%% Measure contrast
[CR,CNR] = measureContrast(b_data_tx,image,-5.5,17.3,4,5,8,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/DAS_cyst_indicated']);

f9 = figure
subplot(211);
bar([CR])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
title('Contrast');
ylabel('CR [dB]');
x_pos = [1 2 3 4 5 6 7];
text(x_pos,double(CR(:)),num2str(round(CR(:)),'%d'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
% subplot(212);
% bar([CNR'])
% set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
% set(gca,'XTickLabel',image.tags)
% ylabel('CNR');
saveas(f9,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/contrast'],'eps2c')
