%clear all;
close all;

filename = [data_path,filesep,'FieldII_STAI_dynamic_range_alt_2.uff'];

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
%%
weights = b_data_weights.get_image('none');
%%
mkdir([ustb_path,filesep,'publications/DynamicRage/figures/simulation/'])

das_img = b_data_das.get_image('none').*weights;  % Compensation weighting
das_img = db(abs(das_img./max(das_img(:))));                 % Normalize on max
f1 = figure(1);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,das_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/DAS'],'eps2c')
axis([-17 2 33 42]);
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/DAS_zoomed'],'eps2c')

cf_img = b_data_cf.get_image('none').*weights;
cf_img = db(abs(cf_img./max(cf_img(:))));                 % Normalize on max
f2 = figure(2);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,cf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f2,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/CF'],'eps2c')
axis([-17 2 33 42]);
saveas(f2,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/CF_zoomed'],'eps2c')

pcf_img = b_data_pcf.get_image('none').*weights;
pcf_img = db(abs(pcf_img./max(pcf_img(:))));
f3 = figure(3);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,pcf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f3,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/PCF'],'eps2c')
axis([-17 2 33 42]);
saveas(f3,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/PCF_zoomed'],'eps2c')

gcf_img = b_data_gcf.get_image('none').*weights;
gcf_img = db(abs(gcf_img./max(gcf_img(:))));
f4 = figure(4);
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,gcf_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f4,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/GCF'],'eps2c')
axis([-17 2 33 42]);
saveas(f4,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/GCF_zoomed'],'eps2c')

%%
mv_img = b_data_mv.get_image('none').*weights;
mv_img = db(abs(mv_img./max(mv_img(:))));
f5 = figure(6);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,mv_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f5,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/MV'],'eps2c')
axis([-17 2 33 42]);
saveas(f5,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/MV_zoomed'],'eps2c')
%%
dmas_img = b_data_dmas.get_image('none').*weights;
dmas_img = db(abs(dmas_img./max(dmas_img(:))));
f7 = figure(7);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,dmas_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f7,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/DMAS'],'eps2c')
axis([-17 2 33 42]);
saveas(f7,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/DMAS_zoomed'],'eps2c')

ebmv_img = b_data_ebmv.get_image('none').*weights;
ebmv_img = db(abs(ebmv_img./max(ebmv_img(:))));
f6 = figure(8);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,ebmv_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f6,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/EBMV'],'eps2c')
axis([-17 2 33 42]);
saveas(f6,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/EBMV_zoomed'],'eps2c')


%%
addpath([ustb_path,'/publications/DynamicRage/functions']);

image.all{1} = das_img;
image.all_intensity{1} = double(b_data_das.get_image('none').*weights);
image.all_intensity{1} = abs(image.all_intensity{1}./max(image.all_intensity{1}(:)));
image.tags{1} = 'DAS';
image.all{5} = mv_img;
image.all_intensity{5} = double(b_data_mv.get_image('none').*weights);
image.all_intensity{5} = abs(image.all_intensity{5}./max(image.all_intensity{5}(:)));
image.tags{5} = 'MV';
image.all{6} = ebmv_img;
image.all_intensity{6} = double(b_data_ebmv.get_image('none').*weights);
image.all_intensity{6} = abs(image.all_intensity{6}./max(image.all_intensity{6}(:)));
image.tags{6} = 'EBMV';
image.all{2} = cf_img;
image.all_intensity{2} = double(b_data_cf.get_image('none').*weights);
image.all_intensity{2} = abs(image.all_intensity{2}./max(image.all_intensity{2}(:)));
image.tags{2} = 'CF';
image.all{7} = pcf_img;
image.all_intensity{7} = double(b_data_pcf.get_image('none').*weights);
image.all_intensity{7} = abs(image.all_intensity{7}./max(image.all_intensity{7}(:)));
image.tags{7} = 'PCF';
image.all{3} = gcf_img;
image.all_intensity{3} = double(b_data_gcf.get_image('none').*weights);
image.all_intensity{3} = abs(image.all_intensity{3}./max(image.all_intensity{3}(:)));
image.tags{3} = 'GCF';
image.all{4} = dmas_img;
image.all_intensity{4} = double(b_data_dmas.get_image('none').*weights);
image.all_intensity{4} = abs(image.all_intensity{4}./max(image.all_intensity{4}(:)));
image.tags{4} = 'DMAS';

%%
[meanLines,x_axis] = getMeanLateralLines(b_data_das,image,47.5,52.5,-20,20);

mask=abs(b_data_tx.scan.x_axis)<10e-3;
theory=-60*(b_data_tx.scan.x_axis(mask)+10e-3)/30e-3-15;

f88 = figure(90);clf;hold all;
subplot(211);
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1});hold on;
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2});
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3});
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4});
plot(b_data_tx.scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color','r');
ylim([-100 0]);
xlim([-20 20])
ylim([-100 0]);
xlim([-20 20])
legend show; grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
title('Lateral gradient');

subplot(212); hold all;
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5});
plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6});
plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7});
plot(b_data_tx.scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color','r');
ylim([-100 0]);
xlim([-20 20])
legend show; grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
%saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral'],'eps2c')

%%
load colorbrewer.mat
colors = colorbrewer.qual.Dark2{8}./255;

[meanLines,x_axis] = getMeanLateralLines(b_data_das,image,47.5,52.5,-12.5,12.5);

mask=abs(b_data_tx.scan.x_axis)<12.5e-3;
theory=-60*(b_data_tx.scan.x_axis(mask)+12.5e-3)/30e-3;

f88 = figure(88);clf;
subplot(211);hold all;
plot(b_data_tx.scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(5,:));
plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(6,:));
plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(7,:));

ylim([-100 0]);
xlim([-12.5 12.5])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral_in_one'],'eps2c')
saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral_in_one'],'png')
%%
%Using "colorbrewer" colors
colors =     [27   158   119; ...
             217    95     2; ...
             117   112   179; ...
             231    41   138; ...
             102   166    30; ...
             230   171     2; ...
             166   118    29; ...
             102   102   102]/255;
%%
[meanLines,x_axis] = getMeanLateralLines(b_data_das,image,47.5,52.5,-12.5,12.5);

mask=abs(b_data_tx.scan.x_axis)<12.5e-3;
theory=-60*(b_data_tx.scan.x_axis(mask)+12.5e-3)/30e-3;

f88 = figure(88);clf;
subplot(211);hold all;
plot(b_data_tx.scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));

ylim([-100 0]);
xlim([-12.5 12.5])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
%title('Lateral gradient');

subplot(212);hold all;
plot(b_data_tx.scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(5,:));
plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(6,:));
plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(7,:));

ylim([-100 0]);
xlim([-12.5 12.5])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral'],'eps2c')
saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral'],'png')
%%

[meanLines_axial,z_axis] = getMeanAxialLines(b_data_das,image,12.5,42.5,11,14);

mask= b_data_tx.scan.z_axis<42.5e-3 & b_data_tx.scan.z_axis>12.5e-3;
theory=-60*(b_data_tx.scan.z_axis(mask)-12.5e-3)/30e-3;

f33 = figure(33);clf; hold all;
subplot(211);hold all;
plot(b_data_tx.scan.z_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(z_axis,meanLines_axial.all{1}-max(meanLines_axial.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(z_axis,meanLines_axial.all{2}-max(meanLines_axial.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(z_axis,meanLines_axial.all{3}-max(meanLines_axial.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(z_axis,meanLines_axial.all{4}-max(meanLines_axial.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));
ylim([-100 0]);
xlim([12.5 42.5]);
legend show; grid on;
ylabel('Normalized amplitude [dB]');
xlabel('z [mm]');
%title('Axial gradient');

subplot(212);hold all;
plot(b_data_tx.scan.z_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(z_axis,meanLines_axial.all{5}-max(meanLines_axial.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(5,:));
plot(z_axis,meanLines_axial.all{6}-max(meanLines_axial.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(6,:));
plot(z_axis,meanLines_axial.all{7}-max(meanLines_axial.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(7,:));

ylim([-100 0]);
xlim([12.5 42.5]);
legend show; grid on;
ylabel('Normalized amplitude [dB]');
xlabel('z [mm]');
%title('Axial gradient');
saveas(f33,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/axial_gradient'],'eps2c')

%%
[meanLines,x_axis] = getMeanLateralLines(b_data_das,image,35,40,b_data_tx.scan.x(1)*10^3,b_data_tx.scan.x(end)*10^3);
theoretical = [-100 ones(1,255)*0 ones(1,255)*-10 -100 ones(1,256)*-100 ones(1,256)*-100 -100 ones(1,255)*0 ones(1,255)*-35 -100];
x_axis_theoretical = linspace(-15,0,1536);

f89 = figure(89);clf;hold all;
subplot(211);hold on;
plot(x_axis_theoretical,theoretical,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));
% ylim([-80 0]);
% xlim([-17.5 1.5]);
% legend('Location','best'); grid on;
% ylabel('Normalized amplitude [dB]');
% xlabel('x [mm]');
% title('Boxes');
% subplot(212); hold all;
%plot(x_axis_theoretical,theoretical,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color','r');
plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(5,:));
plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(6,:));
plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(7,:));
ylim([-80 0]);
xlim([-17.5 1.5]);
legend('Location','best'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
%title('Boxes');
%%
saveas(f89,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/boxes'],'eps2c')


%%

%[handle_1] = plotLateralLine(b_data_tx,image,15,-9,-6,colors)
[handle_2] = plotLateralLine(b_data_tx,image,45,-9,-6,colors)

fwhm_1 = measureFWHM(b_data_tx,image,15,-9,-6)
fwhm_2 = measureFWHM(b_data_tx,image,45,-9,-6)


f9 = figure(9);clf
ax = subplot(2,1,1);
% copyobj(allchild(handle_1),ax)
% ylim([-80 0]);xlim([-9 -6]); 
% legend('show','Location','best');
% ylabel('Amplitude [dB]');xlabel('x [mm]');
% title('PSF at z = 15 mm');
ax = subplot(2,1,1);
copyobj(allchild(handle_2),ax)
ylim([-80 0]);xlim([-9 -6]); 
legend('show','Location','best');
ylabel('Amplitude [dB]');xlabel('x [mm]');
grid on
%title('PSF at z = 45 mm');
saveas(f9,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/PSF_at_z45'],'eps2c')

%%
f10 = figure(10);
subplot(211);
bar([fwhm_1' fwhm_2'])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('FWHM [mm]');
legend('PSF at z = 15 mm','PSF at z = 45 mm')
x_pos = [0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75 7.25];
fwhm_all = [];
for i = 1:length(fwhm_1) 
    fwhm_all = [fwhm_all; fwhm_1(i); fwhm_2(i)];
end
text(x_pos,double(fwhm_all(:)),num2str(fwhm_all(:),'%0.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
grid on;
saveas(f10,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/FWHM'],'eps2c')
%% Measure contrast
[CR,CR_ratio,CNR,C,CNR_picmus,v,v_db] = measureContrast(b_data_tx,image,-7.5,25,3.3,4.5,7.5,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/DAS_cyst_indicated']);

%%

f9 = figure
subplot(411);
bar([CR])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CR [dB]');
x_pos = [1 2 3 4 5 6 7];
text(x_pos,double(CR(:)),num2str(round(CR(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
title('Contrast calculated after log-compression');

subplot(412);
bar([CNR_picmus])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CNR_{picmus}');
x_pos = [1 2 3 4 5 6 7];
text(x_pos,double(CNR_picmus(:)),num2str(round(CNR_picmus(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)

subplot(413);
bar([C])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('C');
x_pos = [1 2 3 4 5 6 7];
text(x_pos,double(C(:)),num2str((C(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)

title('Contrast calculated before log-compression');
subplot(414);
bar([CNR])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CNR');
x_pos = [1 2 3 4 5 6 7];
text(x_pos,double(CNR(:)),num2str((CNR(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)

%%

%saveas(f9,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/contrast_all'],'eps2c')
%%
f999 = figure(999);
subplot(211);
bar(v_db)
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('variance');
x_pos = [1 2 3 4 5 6 7];
text(x_pos,double(v_db(:)),num2str(round(v_db(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
title('Variance of speckle after log-compression');

subplot(212);
bar(v)
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('variance');
x_pos = [1 2 3 4 5 6 7];
text(x_pos,double(v(:)),num2str((v(:)),'%.6f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
title('Variance of speckle before log-compression');


saveas(f999,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/variance'],'eps2c')
%%

figure;
subplot(421)
imagesc(image.all{1})
caxis([-60 0])
colorbar
title(image.tags{1})
subplot(422)
imagesc(image.all_intensity{1})
colorbar
title(image.tags{1})
subplot(423)
imagesc(image.all{2})
colorbar
caxis([-60 0])
title(image.tags{2});
subplot(424)
imagesc(image.all_intensity{2})
colorbar
title(image.tags{2});
subplot(425)
imagesc(image.all{6})
colorbar
caxis([-60 0])
title(image.tags{6});
subplot(426)
imagesc(image.all_intensity{6})
colorbar
title(image.tags{6});
subplot(427)
imagesc(image.all{5})
colorbar
caxis([-60 0])
title(image.tags{5});
subplot(428)
imagesc(image.all_intensity{5})
colorbar
title(image.tags{5});
