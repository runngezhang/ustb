% addpath make_simulation
% 
% [point_position, point_amplitudes] = simulatedPhantomDynamicRange_2(650);
% 
% %%
% figure(10);
% plot3(point_position(:,1)*1e3,point_position(:,3)*1e3,20*log10(point_amplitudes),'b.'); axis equal;
% zlim([-60 5]);
% xlabel('X [mm]');ylabel('Z [mm]');zlabel('Amplitude [dB]');
% 
% %%
% addpath([ustb_path,filesep,'publications/DynamicRage/functions'])
% addpath([ustb_path,filesep,'publications/DynamicRage/functions/tightfig'])
% f1 = figure(11);
% scatter3(point_position(:,1)*1e3,point_position(:,3)*1e3,20*log10(point_amplitudes),13,20*log10(point_amplitudes),'filled'); axis equal;
% colormap gray;
% 
% az = 0;
% el = 90;
% view(az, el);
% %set(gca, 'Xdir', 'reverse')
% set(gca, 'Ydir', 'reverse')
% xlabel('x [mm]');ylabel('z [mm]');
% 
% set(gca,'Color','k')
% ylim([10 55])
% %xlim([-20 20])
% caxis([-60 0])
% 
% axi = gca;
% 
% xc_nonecho = -7.5;
% zc_nonecho = 25;
% r_nonecho = 3.3;
% r_speckle_inner = 4.5;
% r_speckle_outer = 7.5;
% viscircles(axi,[xc_nonecho,zc_nonecho],r_nonecho,'EdgeColor','r','EnhanceVisibility',0);
% viscircles(axi,[xc_nonecho,zc_nonecho],r_speckle_inner,'EdgeColor','b','EnhanceVisibility',0);
% viscircles(axi,[xc_nonecho,zc_nonecho],r_speckle_outer,'EdgeColor','b','EnhanceVisibility',0);
% 
% set(gca,'FontSize',14)
% axis image;
% f1 = tightfig(f1);
% f1.InvertHardcopy = 'off';
% saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/GT_alt_1'],'eps2c')
% savefig(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/GT_alt_1'])
% %%
% caxis([-90 0])
% saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/GT_alt_2'],'eps2c')
% 
% h = colorbar;
% ylabel(h,'Amplitude [dB]');
% saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/GT_alt_3'],'eps2c')
% caxis([-60 0])
% saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/GT_alt_4'],'eps2c')
% 
% %%
% axis([-17 2 33 42]);
% caxis([-60 0])
% saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/GT_zoomed'],'eps2c')
% 
% 
% %%
% 
% Z = [point_position(:,1)*1e3 point_position(:,3)*1e3 (point_amplitudes)];
% 
% %%
% figure;
% 
% mesh(Z); axis equal;

scan=uff.linear_scan('x_axis',linspace(-20e-3,20e-3,1024).', 'z_axis', linspace(8e-3,55e-3,2048).');

b_data_t = uff.beamformed_data();
b_data_t.scan = scan;

%%

mask_x_speckle_with_cyst = logical(b_data_t.scan.x_axis > -15e-3) & logical(b_data_t.scan.x_axis < 0);
mask_z_speckle_with_cyst = logical(b_data_t.scan.z_axis > 17.5e-3) & logical(b_data_t.scan.z_axis < 32.5e-3);

mask_x_gradient_lateral= logical(b_data_t.scan.x_axis > -20e-3) & logical(b_data_t.scan.x_axis < 20e-3);
mask_z_gradient_lateral = logical(b_data_t.scan.z_axis > 47.5e-3) & logical(b_data_t.scan.z_axis < 52.5e-3);

mask_x_gradient_axial = logical(b_data_t.scan.x_axis > 10e-3) & logical(b_data_t.scan.x_axis < 15e-3);
mask_z_gradient_axial = logical(b_data_t.scan.z_axis > 10e-3) & logical(b_data_t.scan.z_axis < 45e-3);

db_fall_per_mm = 2;
theory_lateral_gradient_amp= -(40*db_fall_per_mm)*(b_data_t.scan.x_axis(mask_x_gradient_lateral)+20e-3)/40e-3;
z_pos = (b_data_t.scan.z_axis(mask_z_gradient_axial)-min(b_data_t.scan.z_axis(mask_z_gradient_axial)));
z_pos = z_pos./max(z_pos);
theory_axial_gradient_amp= -(35*db_fall_per_mm)*z_pos;

mask_x_box_15_12 = logical(b_data_t.scan.x_axis > -15e-3) & logical(b_data_t.scan.x_axis < -12.5e-3);
mask_x_box_12_10 = logical(b_data_t.scan.x_axis > -12.5e-3) & logical(b_data_t.scan.x_axis < -10e-3);
mask_x_box_5_25 = logical(b_data_t.scan.x_axis > -5e-3) & logical(b_data_t.scan.x_axis < -2.5e-3);
mask_x_box_25_0 = logical(b_data_t.scan.x_axis > -2.5e-3) & logical(b_data_t.scan.x_axis < -0e-3);
mask_z_box = logical(b_data_t.scan.z_axis > 35e-3) & logical(b_data_t.scan.z_axis < 40e-3);


% Create masks to mask out the ROI of the cyst and the background.
for p = 1:length(b_data_t.scan.z_axis)
    positions(p,:,1) = b_data_t.scan.x_axis;
end
for p = 1:length(b_data_t.scan.x_axis)
    positions(:,p,2) = b_data_t.scan.z_axis;
end

% Coordinates for cyst
xc_cyst = -7.5*10^-3;
zc_cyst = 25*10^-3;
r_cyst = 4*10^-3;
points = ((positions(:,:,1)-xc_cyst).^2) + (positions(:,:,2)-zc_cyst).^2;
idx_cyst = (points < (r_cyst)^2);                     %ROI inside cyst

% Coordinates for PSF
xc_psf = -7.5*10^-3;
zc_psf_1 = 15*10^-3;
zc_psf_2 = 45*10^-3;
r_psf = 0.25*10^-3;

points = ((positions(:,:,1)-xc_psf).^2) + (positions(:,:,2)-zc_psf_1).^2;
idx_psf_1 = (points < (r_psf)^2);                     %ROI inside cyst
points = ((positions(:,:,1)-xc_psf).^2) + (positions(:,:,2)-zc_psf_2).^2;
idx_psf_2 = (points < (r_psf)^2);                     %ROI inside cyst

theoretical_image = ones(b_data_t.scan.N_z_axis,b_data_t.scan.N_x_axis)*-1000;
theoretical_image(mask_z_speckle_with_cyst,mask_x_speckle_with_cyst) = -10;
theoretical_image(mask_z_gradient_lateral,mask_x_gradient_lateral) = repmat(theory_lateral_gradient_amp',[size(b_data_t.scan.z_axis(mask_z_gradient_lateral)),1]);
theoretical_image(mask_z_gradient_axial,mask_x_gradient_axial) = repmat(theory_axial_gradient_amp,[1,size(b_data_t.scan.x_axis(mask_x_gradient_axial))]);
theoretical_image(mask_z_box,mask_x_box_15_12) = 0;
theoretical_image(mask_z_box,mask_x_box_12_10) = -15;
theoretical_image(mask_z_box,mask_x_box_5_25) = 0;
theoretical_image(mask_z_box,mask_x_box_25_0) = -35;
theoretical_image(idx_cyst(:)) = -1000;
theoretical_image(idx_psf_1(:)) = 0;
theoretical_image(idx_psf_2(:)) = 0;


imagesc(theoretical_image)
caxis([-60 0])
colorbar

% Make it signal strength not dB
theoretical_image_non_db = 10.^(theoretical_image/20);


b_data_t.data = theoretical_image_non_db(:);
b_data_t.plot([],[],60)

%%

theory_img = b_data_t.get_image('none');  % Compensation weighting
theory_img = db(abs(theory_img./max(theory_img(:))));                 % Normalize on max
f1 = figure(1);clf;
imagesc(b_data_t.scan.x_axis*1000,b_data_t.scan.z_axis*1000,theory_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/theoretical'],'eps2c')
axis([-17 2 33 42]);
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/simulation/theoretical_zoomed'],'eps2c')
