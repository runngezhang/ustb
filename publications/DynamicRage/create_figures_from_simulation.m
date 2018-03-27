clear all;
close all;
%%
filename = [data_path,filesep,'FieldII_STAI_dynamic_range_alt_2.uff'];

b_data_tx = uff.beamformed_data();
b_data_das = uff.beamformed_data();
b_data_cf = uff.beamformed_data();
b_data_pcf = uff.beamformed_data();
b_data_gcf = uff.beamformed_data();
b_data_mv = uff.beamformed_data();
b_data_ebmv = uff.beamformed_data();
b_data_dmas = uff.beamformed_data();
b_data_glt = uff.beamformed_data();
b_data_weights = uff.beamformed_data();

b_data_tx.read(filename,'/b_data_tx');
b_data_das.read(filename,'/b_data_das');
b_data_cf.read(filename,'/b_data_cf');
b_data_pcf.read(filename,'/b_data_pcf');
b_data_gcf.read(filename,'/b_data_gcf');
b_data_mv.read(filename,'/b_data_mv');
b_data_ebmv.read(filename,'/b_data_ebmv');
b_data_dmas.read(filename,'/b_data_dmas');
b_data_glt.read(filename,'/b_data_glt');
b_data_weights.read(filename,'/b_data_weights');

weights = b_data_weights.get_image('none');
%%
mkdir([ustb_path,filesep,'publications/DynamicRage/figures/simulation/'])

das_img = b_data_das.get_image('none').*weights;  % Compensation weighting
das_img = db(abs(das_img./max(das_img(:))));      % Normalize on max
f1 = figure(1);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,das_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/DAS'],'eps2c')
axis([-17 2 33 42]);
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/DAS_zoomed'],'eps2c')


% power_das_img = (b_data_das.get_image('none').*weights).^1.4;  % Compensation weighting
% power_das_img = db(abs(power_das_img./max(power_das_img(:))));      % Normalize on max
% f10 = figure(10);clf;
% imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,power_das_img);
% colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
% set(gca,'FontSize',14)
% saveas(f10,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/PDAS'],'eps2c')
% axis([-17 2 33 42]);
% saveas(f10,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/PDAS_zoomed'],'eps2c')
% 
% %% gamma function we want to approximate
% a=0.00015;
% b=0;
% c=0.93;
% 
% a = 0;
% b = 0.01;
% c = 0.8;
% 
% 
% % linear space
% x=logspace(-60/20,0,100);
% 
% % dB space
% x_dB=20*log10(x);
% %x_dB_compressed=-a*x_dB.^2+b*x_dB;
% x_dB_compressed=a*x_dB.^3-b*x_dB.^2+c*x_dB;
% x_dB_compressed_ref = -0.01*x_dB.^2+0.8*x_dB;
% % find the cublic spline that approximate the compressed values
% x_compressed=10.^(x_dB_compressed/20);
% x_compressed_ref=10.^(x_dB_compressed_ref/20);
% gamma = fit(x.',x_compressed.','cubicspline');
% gamma_ref = fit(x.',x_compressed_ref.','cubicspline');
% 
% signal = abs(b_data_das.get_image('none').*weights);
% signal = signal./max(signal(:));
% gamma_img = gamma(signal);
% gamma_img_signal = reshape(gamma_img,2048,1024);
% gamma_img = db(gamma_img_signal./max(gamma_img_signal(:)));
% 
% figure(8888);clf;
% subplot(1,2,1);
% plot(x,x,'k'); hold on; grid on; axis equal tight;
% plot(x,x_compressed,'b'); hold on;
% plot(x,gamma(x),'r:');
% %plot(x,gamma_ref(x),'b--');
% title('Linear space');
% 
% subplot(1,2,2);hold all;
% plot(x_dB,x_dB,'k'); hold on; grid on; axis equal tight;
% plot(x_dB,x_dB_compressed,'b'); hold on;
% plot(x_dB,20*log10(gamma(x)),'--'); hold on;
% %plot(x_dB,20*log10(gamma_ref(x)),'k--'); hold on;
% title('Log space');
%%
glt_img = b_data_glt.get_image('none').*weights;  % Compensation weighting
glt_img = db(abs(glt_img./max(glt_img(:))));      % Normalize on max

f100 = figure(100);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,glt_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f100,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/GLT'],'eps2c')
axis([-17 2 33 42]);
saveas(f100,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/GLT_zoomed'],'eps2c')

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

mv_img = b_data_mv.get_image('none').*weights;
mv_img = db(abs(mv_img./max(mv_img(:))));
f5 = figure(6);clf;
imagesc(b_data_das.scan.x_axis*1000,b_data_das.scan.z_axis*1000,mv_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14)
saveas(f5,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/MV'],'eps2c')
axis([-17 2 33 42]);
saveas(f5,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/MV_zoomed'],'eps2c')

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
addpath([ustb_path,filesep,'publications/DynamicRage/functions/tightfig'])

image.all{1} = das_img;
image.all_signal{1} = double(b_data_das.get_image('none').*weights);
image.all_signal{1} = (image.all_signal{1}./max(image.all_signal{1}(:)));
image.tags{1} = 'DAS';
image.all{5} = mv_img;
image.all_signal{5} = double(b_data_mv.get_image('none').*weights);
image.all_signal{5} = (image.all_signal{5}./max(image.all_signal{5}(:)));
image.tags{5} = 'MV';
image.all{6} = ebmv_img;
image.all_signal{6} = double(b_data_ebmv.get_image('none').*weights);
image.all_signal{6} = (image.all_signal{6}./max(image.all_signal{6}(:)));
image.tags{6} = 'EBMV';
image.all{2} = cf_img;
image.all_signal{2} = double(b_data_cf.get_image('none').*weights);
image.all_signal{2} = (image.all_signal{2}./max(image.all_signal{2}(:)));
image.tags{2} = 'CF';
image.all{7} = pcf_img;
image.all{7}(695,320) = image.all{7}(694,320);
image.all_signal{7} = double(b_data_pcf.get_image('none').*weights);
image.all_signal{7} = (image.all_signal{7}./max(image.all_signal{7}(:)));
image.all_signal{7}(695,320) = image.all_signal{7}(694,320);
image.tags{7} = 'PCF';
image.all{3} = gcf_img;
image.all_signal{3} = double(b_data_gcf.get_image('none').*weights);
image.all_signal{3} = (image.all_signal{3}./max(image.all_signal{3}(:)));
image.tags{3} = 'GCF';
image.all{4} = dmas_img;
image.all_signal{4} = double(b_data_dmas.get_image('none').*weights);
image.all_signal{4} = (image.all_signal{4}./max(image.all_signal{4}(:)));
image.tags{4} = 'DMAS';
image.all{8} = glt_img;
image.all_signal{8} = double(b_data_glt.get_image('none').*weights);
image.all_signal{8} = (image.all_signal{8}./max(image.all_signal{8}(:)));
image.tags{8} = 'GLT';
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
% [meanLines,x_axis] = getMeanLateralLines(b_data_das,image,47.5,52.5,-12.5,12.5);
% 
% mask=abs(b_data_tx.scan.x_axis)<12.5e-3;
% theory=-60*(b_data_tx.scan.x_axis(mask)+12.5e-3)/30e-3;
% 
% f88 = figure(88);clf;
% subplot(211);hold all;
% plot(b_data_tx.scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
% plot(x_axis,meanLines.all{1}-max(meanLines.all{1}),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
% plot(x_axis,meanLines.all{2}-max(meanLines.all{2}),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
% plot(x_axis,meanLines.all{3}-max(meanLines.all{3}),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
% plot(x_axis,meanLines.all{4}-max(meanLines.all{4}),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));
% plot(x_axis,meanLines.all{5}-max(meanLines.all{5}),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(5,:));
% plot(x_axis,meanLines.all{6}-max(meanLines.all{6}),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(6,:));
% plot(x_axis,meanLines.all{7}-max(meanLines.all{7}),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(7,:));
% 
% ylim([-100 0]);
% xlim([-12.5 12.5])
% legend('Location','sw'); grid on;
% ylabel('Normalized amplitude [dB]');
% xlabel('x [mm]');
% saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral_in_one'],'eps2c')
% saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral_in_one'],'png')
%%
%Using "colorbrewer" colors
colors =     [27   158   119; ...
             217    95     2; ...
             117   112   179; ...
             231    41   138; ...
             102   166    30; ...
             230   171     2; ...
             166   118    29; ...
             102   102   102; ...
             0      0      0]/255;
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
plot(x_axis,meanLines.all{8}-max(meanLines.all{8}),'Linewidth',2,'DisplayName',image.tags{8},'Color',colors(8,:));


ylim([-100 0]);
xlim([-12.5 12.5])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral'],'eps2c')
saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral'],'png')
savefig(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral'])
%%
image_experimental = load('Experimental.mat','image');
image_experimental = image_experimental.image;
b_data_das_exp = load('Experimental.mat','b_data_das'); 
b_data_das_exp = b_data_das_exp.b_data_das;

[meanLinesLateralExp,~] = getMeanLateralLines(b_data_das_exp,image_experimental,39,48.5,-12.5,12.5);


mask_lateral_exp=abs(b_data_das_exp.scan.x_axis)<12.5e-3;
theory_lateral_exp=-40*(b_data_das_exp.scan.x_axis(mask_lateral_exp)+12.5e-3)/25e-3;



[meanLinesLateral,x_axis] = getMeanLateralLines(b_data_das,image,47.5,52.5,-10,10);
mask_lateral=abs(b_data_tx.scan.x_axis)<10e-3;
theory_lateral=-40*(b_data_tx.scan.x_axis(mask_lateral)+10e-3)/20e-3;

[meanLines_axial,z_axis] = getMeanAxialLines(b_data_das,image,20,40,11,14);
mask_axial= b_data_tx.scan.z_axis<40e-3 & b_data_tx.scan.z_axis>20e-3;
theory_axial=-40*(b_data_tx.scan.z_axis(mask_axial)-20e-3)/20e-3;

for i = 1:length(image.all)
    f88 = figure(8888+i);clf;hold all;
    plot(theory_lateral,meanLinesLateral.all{i}-max(meanLinesLateral.all{i}),'LineStyle','-.','Linewidth',2,'DisplayName','Sim. lateral');
    plot(theory_axial,meanLines_axial.all{i}-max(meanLines_axial.all{i}),'Linestyle','--','Linewidth',2,'DisplayName','Sim. axial');
    plot(theory_lateral_exp,meanLinesLateralExp.all{i}-max(meanLinesLateralExp.all{i}),'Linestyle',':','Linewidth',2,'DisplayName','Exp. lateral');
    plot(theory_lateral,theory_lateral,'LineStyle','-','Linewidth',2,'Color','k','DisplayName','Theory');
    set(gca, 'XDir','reverse');%title(image.tags{i});
    axis image
    legend('show','Location','sw');
    ylim([-80 0])
    xlim([-40 0])
    grid on
    xlabel('Theoretical [dB]');
    ylabel('Output [dB]');
    saveas(f88,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient/',image.tags{i}],'eps2c')
end



%%
f900 = figure(900);clf;
subplot(211);hold all;
plot(b_data_tx.scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(x_axis,db(meanLines.all_signal{1}./max(meanLines.all_signal{1})),'Linewidth',2,'DisplayName',image.tags{1},'Color',colors(1,:));
plot(x_axis,db(meanLines.all_signal{2}./max(meanLines.all_signal{2})),'Linewidth',2,'DisplayName',image.tags{2},'Color',colors(2,:));
plot(x_axis,db(meanLines.all_signal{3}./max(meanLines.all_signal{3})),'Linewidth',2,'DisplayName',image.tags{3},'Color',colors(3,:));
plot(x_axis,db(meanLines.all_signal{4}./max(meanLines.all_signal{4})),'Linewidth',2,'DisplayName',image.tags{4},'Color',colors(4,:));

ylim([-100 0]);
xlim([-12.5 12.5])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
%title('Lateral gradient');
%%
subplot(212);hold all;
plot(b_data_tx.scan.x_axis(mask)*10^3,theory,'LineStyle','-.','Linewidth',2,'DisplayName','Theoretical','Color',[0 0 0]);
plot(x_axis,db(meanLines.all_signal{5}./max(meanLines.all_signal{5})),'Linewidth',2,'DisplayName',image.tags{5},'Color',colors(5,:));
plot(x_axis,db(meanLines.all_signal{6}./max(meanLines.all_signal{6})),'Linewidth',2,'DisplayName',image.tags{6},'Color',colors(6,:));
plot(x_axis,db(meanLines.all_signal{7}./max(meanLines.all_signal{7})),'Linewidth',2,'DisplayName',image.tags{7},'Color',colors(7,:));

ylim([-100 0]);
xlim([-12.5 12.5])
legend('Location','sw'); grid on;
ylabel('Normalized amplitude [dB]');
xlabel('x [mm]');
saveas(f900,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral_mean_signal'],'eps2c')
saveas(f900,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/gradient_lateral_mean_signal'],'png')


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
plot(x_axis,meanLines.all{8}-max(meanLines.all{8}),'Linewidth',2,'DisplayName',image.tags{8},'Color',colors(8,:));
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
%%
saveas(f9,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/PSF_at_z45'],'eps2c')

%%
f10 = figure(10);clf;
subplot(211);
bar([fwhm_1' fwhm_2'])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('FWHM [mm]');
legend('Location','best','PSF at z = 15 mm','PSF at z = 45 mm')
x_pos = [0.75 1.25 1.75 2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75 7.25 7.75 8.25 8.75 9.25];
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
[CR_signal, CR_signal_dagger, CR_image, CNR_signal, CNR_image] = measureContrast(b_data_tx,image,-7.5,25,3.3,4.5,7.5,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/DAS_cyst_indicated']);

%%

f9 = figure
subplot(211);
bar([10*log10(CR_signal)])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CR_{signal} [dB]');
x_pos = [1 2 3 4 5 6 7 8];
text(x_pos,double(10*log10(CR_signal(:))),num2str(round(10*log10(CR_signal(:))),'%d'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
set(gca,'Ydir','reverse')

subplot(212);
bar([CR_image])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CR_{image} [dB]');
x_pos = [1 2 3 4 5 6 7 8];
text(x_pos,double(CR_image(:)),num2str(round(CR_image(:)),'%d'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)

saveas(f9,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/CR'],'eps2c')

saveas(f9,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/CR'],'png')
%%

f77 = figure(77);
subplot(211);
bar([CNR_signal])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CNR_{signal}');
x_pos = [1 2 3 4 5 6 7 8];
text(x_pos,double(CNR_signal(:)),num2str((CNR_signal(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)
ylim([0 1.1])

subplot(212);
bar([CNR_image])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CNR_{image}');
x_pos = [1 2 3 4 5 6 7 8];
text(x_pos,double(CNR_image(:)),num2str((CNR_image(:)),'%.2f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)

saveas(f77,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/CNR'],'eps2c')
saveas(f77,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/CNR'],'png')
%%
f99 = figure(99);
subplot(411);
bar([CR_signal_dagger])
set(gca,'XTick',linspace(1,length(image.tags),length(image.tags)))
set(gca,'XTickLabel',image.tags)
ylabel('CR power');
x_pos = [1 2 3 4 5 6 7 8];
text(x_pos,double(CR_signal_dagger(:)-0.2),num2str((CR_signal_dagger(:)),'%.6f'),...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',12)

ylim([-1.3 0])