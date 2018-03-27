addpath make_simulation

[point_position, point_amplitudes] = mimicing_experimentalPhantomDynamicRange(650);

%%
figure(10);
plot3(point_position(:,1)*1e3,point_position(:,3)*1e3,20*log10(point_amplitudes),'b.'); axis equal;
zlim([-60 5]);
xlabel('X [mm]');ylabel('Z [mm]');zlabel('Amplitude [dB]');

%%
addpath([ustb_path,filesep,'publications/DynamicRage/functions'])
addpath([ustb_path,filesep,'publications/DynamicRage/functions/tightfig'])
f1 = figure(11);
scatter3(point_position(:,1)*1e3,point_position(:,3)*1e3,20*log10(point_amplitudes),13,20*log10(point_amplitudes),'filled'); axis equal;
colormap gray;

az = 0;
el = 90;
view(az, el);
set(gca, 'Ydir', 'reverse')
xlabel('x [mm]');ylabel('z [mm]');

set(gca,'Color','k')
ylim([7.5 52])
xlim([-20 20])
caxis([-60 0])

[x, y, A] = getimage(f1);

%%
axi = gca;

xc_nonecho = -5.5;
zc_nonecho = 17.5;
r_nonecho = 4;
r_speckle_inner = 5;
r_speckle_outer = 8;
viscircles(axi,[xc_nonecho,zc_nonecho],r_nonecho,'EdgeColor','r','EnhanceVisibility',0);
viscircles(axi,[xc_nonecho,zc_nonecho],r_speckle_inner,'EdgeColor','b','EnhanceVisibility',0);
viscircles(axi,[xc_nonecho,zc_nonecho],r_speckle_outer,'EdgeColor','b','EnhanceVisibility',0);



%%
set(gca,'FontSize',14)
axis image;
f1 = tightfig(f1);
ylim([7.5 52])
xlim([-20 20])
f1.InvertHardcopy = 'off';
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/GT_alt_1'],'eps2c')
savefig(f1,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/GT_alt_1'])


%%

scan=uff.linear_scan('x_axis',linspace(-20e-3,20e-3,512*2).','z_axis',linspace(8e-3,52e-3,2048).');

b_data_t = uff.beamformed_data();
b_data_t.scan = scan;

%%

mask_x_speckle = logical(b_data_t.scan.x_axis > -15.5e-3) & logical(b_data_t.scan.x_axis < 4.5e-3);
mask_z_speckle = logical(b_data_t.scan.z_axis > 10e-3) & logical(b_data_t.scan.z_axis < 27e-3);

mask_x_gradient= logical(b_data_t.scan.x_axis > -20e-3) & logical(b_data_t.scan.x_axis < 20e-3);
mask_z_gradient = logical(b_data_t.scan.z_axis > 40e-3) & logical(b_data_t.scan.z_axis < 50e-3);

db_fall_per_mm = 1.66;
theory_gradient_amp= -(40*db_fall_per_mm)*(b_data_t.scan.x_axis(mask_x_gradient)+20e-3)/40e-3;

% Create masks to mask out the ROI of the cyst and the background.
for p = 1:length(b_data_t.scan.z_axis)
    positions(p,:,1) = b_data_t.scan.x_axis;
end
for p = 1:length(b_data_t.scan.x_axis)
    positions(:,p,2) = b_data_t.scan.z_axis;
end

% Coordinates for cyst
xc_cyst = -5.5*10^-3;
zc_cyst = 18.75*10^-3;
r_cyst = 8.5/2*10^-3;
points = ((positions(:,:,1)-xc_cyst).^2) + (positions(:,:,2)-zc_cyst).^2;
idx_cyst = (points < (r_cyst)^2);                     %ROI inside cyst

% Coordinates for PSF
xc_psf = 10*10^-3;
zc_psf_1 = 10*10^-3;
zc_psf_2 = 20*10^-3;
zc_psf_3 = 30*10^-3;
r_psf = 0.25*10^-3;

points = ((positions(:,:,1)-xc_psf).^2) + (positions(:,:,2)-zc_psf_1).^2;
idx_psf_1 = (points < (r_psf)^2);                     %ROI inside cyst
points = ((positions(:,:,1)-xc_psf).^2) + (positions(:,:,2)-zc_psf_2).^2;
idx_psf_2 = (points < (r_psf)^2);                     %ROI inside cyst
points = ((positions(:,:,1)-xc_psf).^2) + (positions(:,:,2)-zc_psf_3).^2;
idx_psf_3 = (points < (r_psf)^2);                     %ROI inside cyst

theoretical_image = ones(b_data_t.scan.N_z_axis,b_data_t.scan.N_x_axis)*-1000;
theoretical_image(mask_z_speckle,mask_x_speckle) = -10;
theoretical_image(mask_z_gradient,mask_x_gradient) = repmat(theory_gradient_amp',[size(b_data_t.scan.z_axis(mask_z_gradient)),1]);
theoretical_image(idx_cyst(:)) = -1000;
theoretical_image(idx_psf_1(:)) = 0;
theoretical_image(idx_psf_2(:)) = 0;
theoretical_image(idx_psf_3(:)) = 0;

imagesc(theoretical_image)
caxis([-60 0])
colorbar

% Make it signal strength not dB
theoretical_image_non_db = 10.^(theoretical_image/20);

%%
b_data_t.data = theoretical_image_non_db(:);
b_data_t.plot([],[],60)


theory_img = b_data_t.get_image('none');  % Compensation weighting
theory_img = db(abs(theory_img./max(theory_img(:))));                 % Normalize on max
f1 = figure(1);clf;
imagesc(b_data_t.scan.x_axis*1000,b_data_t.scan.z_axis*1000,theory_img);
colormap gray;caxis([-60 0]);axis image;xlabel('x [mm]');ylabel('z [mm]');
set(gca,'FontSize',14);
saveas(f1,[ustb_path,filesep,'publications/DynamicRage/figures/experimental/theoretical'],'eps2c')
