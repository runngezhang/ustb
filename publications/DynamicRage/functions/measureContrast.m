function [CR_signal, CR_signal_dagger, CR_image, CNR_signal, CNR_image] = measureContrast(sta_image,image,xc_nonecho,zc_nonecho,r_nonecho,r_speckle_inner,r_speckle_outer,f_filename)
%PLOTLATERALLINE Plot lateral line from all images saved in image struct
%   
%% Non echo Cyst contrast

xc_speckle = xc_nonecho;
zc_speckle = zc_nonecho;

% Create masks to mask out the ROI of the cyst and the background.
for p = 1:length(sta_image.scan.z_axis)
    positions(p,:,1) = sta_image.scan.x_axis;
end

for p = 1:length(sta_image.scan.x_axis)
    positions(:,p,2) = sta_image.scan.z_axis;
end


points = ((positions(:,:,1)-xc_nonecho*10^-3).^2) + (positions(:,:,2)-zc_nonecho*10^-3).^2;
idx_cyst = (points < (r_nonecho*10^-3)^2);                     %ROI inside cyst
idx_speckle_outer =  (((positions(:,:,1)-xc_speckle*10^-3).^2) + (positions(:,:,2)-zc_speckle*10^-3).^2 < (r_speckle_outer*10^-3)^2); %ROI speckle
idx_speckle_inner =  (((positions(:,:,1)-xc_speckle*10^-3).^2) + (positions(:,:,2)-zc_speckle*10^-3).^2 < (r_speckle_inner*10^-3)^2); %ROI speckle
idx_speckle_outer(idx_speckle_inner) = 0;% & idx_speckle_inner;
idx_speckle = idx_speckle_outer;

%%

f1 = figure(1111);clf;
%set(f1,'Position',[100,100,300,300]);
%subplot(1,length(image.all),i)
imagesc(sta_image.scan.x_axis*1e3,sta_image.scan.z_axis*1e3,image.all{1});
colormap gray; caxis([-60 0]); axis image; xlabel('x [mm]');ylabel('z [mm]');
%colorbar
axi = gca;
set(gca,'FontSize',14)
viscircles(axi,[xc_nonecho,zc_nonecho],r_nonecho,'EdgeColor','r','EnhanceVisibility',0);
viscircles(axi,[xc_nonecho,zc_nonecho],r_speckle_inner,'EdgeColor','y','EnhanceVisibility',0);
viscircles(axi,[xc_nonecho,zc_nonecho],r_speckle_outer,'EdgeColor','y','EnhanceVisibility',0);
%title(image.tags{1});
%viscircles(axi,[-xc_nonecho,zc_nonecho],r_nonecho,'EdgeColor','r');
%viscircles(axi,[-xc_nonecho,zc_nonecho],r_speckle_inner,'EdgeColor','y');
%viscircles(axi,[-xc_nonecho,zc_nonecho],r_speckle_outer,'EdgeColor','y');

if nargin == 8
    set(gca,'FontSize',15)
    saveas(f1,f_filename,'eps2c');
end

%%

%%
for i = 1:length(image.all)
    u_ROI_signal = mean( abs(image.all_signal{i}(idx_cyst)).^2 );
    u_B_signal = mean( abs(image.all_signal{i}(idx_speckle)).^2 ); 
    
    sigma_ROI_signal = std(abs(image.all_signal{i}(idx_cyst)).^2);
    sigma_B_signal  = std(abs(image.all_signal{i}(idx_speckle)).^2);
    
    u_ROI_image = mean( image.all{i}(idx_cyst) );
    u_B_image = mean( image.all{i}(idx_speckle) ); 
    
    sigma_ROI_image = std( image.all{i}(idx_cyst) );
    sigma_B_image  = std( image.all{i}(idx_speckle) );
    
    %Current definition on the manuscript
    CR_signal(i) =  u_ROI_signal / u_B_signal;
    CR_signal_dagger(i) = (u_ROI_signal - u_B_signal) / sqrt(u_ROI_signal^2+u_B_signal^2);

    CR_image(i) = abs(u_ROI_image - u_B_image);
    
    CNR_signal(i) = abs(u_ROI_signal - u_B_signal) / sqrt((sigma_ROI_signal^2 + sigma_B_signal^2));
    CNR_image(i) = abs(u_ROI_image - u_B_image) / sqrt((sigma_ROI_image^2 + sigma_B_image^2));
end

end

