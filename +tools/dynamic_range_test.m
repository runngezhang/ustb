function [filtered_p,h,w] = dynamic_range_test(channel_data,b_data)


if strcmp(channel_data.name,'Simulated dynamic range phantom. Created with Field II. See the reference for details')
    z_start = 47.5;
    z_stop = 52.5;
    x_start = -15;
    x_stop = 15;
    
    gradient = -2% db/mm
 %   theory_lateral = -40*(b_data.scan.x_axis(mask_lateral)+10e-3)/20e-3;
elseif strcmp(channel_data.name,'Experimental dynamic range phantom. Created with Field II. See the reference for details')
        
else
    error('The dynamic range test is only defined for the simulated dynamic range phantom, and the experimental one. See the example.');
end


% Mask out the top gradient part of the image
mask=reshape(b_data.scan.z>z_start*10^-3&b_data.scan.z<z_stop*10^-3,[length(b_data.scan.z_axis) length(b_data.scan.x_axis)])...
     &reshape(b_data.scan.x>x_start*10^-3&b_data.scan.x<x_stop*10^-3,[length(b_data.scan.z_axis) length(b_data.scan.x_axis)]);

 %%
figure(1211);
imagesc(b_data.scan.x_axis*1e3,b_data.scan.z_axis*1e3,mask.*b_data.get_image());

%%

% Find the x min and max index in the image from the gradient mask, top
z_min_top = rem(min(find(mask==1)),length(b_data.scan.z_axis));
z_max_top = rem(max(find(mask==1)),length(b_data.scan.z_axis));

temp = find(mask(z_min_top,:)==1);
x_min = temp(1);
x_max = temp(end);

img = b_data.get_image();
lateral_line = mean(img(z_min_top:z_max_top,x_min:x_max),1);
lateral_line = lateral_line-max(lateral_line(:)); %Normalize to have the gradient start at 0

regresion_coeff = polyfit(b_data.scan.x_axis(x_min:x_max)*10^3,lateral_line',1);
regression = polyval(regresion_coeff,b_data.scan.x_axis(x_min:x_max)*10^3);

theory_lateral = (gradient*(b_data.scan.x_axis(x_min:x_max)-x_start*10^-3))*10^3;

figure
hold all;
plot(theory_lateral,theory_lateral);
plot(theory_lateral,lateral_line);
plot(theory_lateral,regression);
set(gca, 'XDir','reverse');
title(sprintf('Theoretical gradient: %.4f dB/mm, estimated: %.4f dB/mm',gradient,regresion_coeff(1)));
end

