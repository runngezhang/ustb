function [image_corrected, correcting_coeff] = calibrateImage(sta_image,image,z_start,z_stop,x_start,x_stop,order)
%CALIBRATEIMAGE Summary of this function goes here
%   Detailed explanation goes here

[meanLines,x_axis] = getMeanLateralLines(sta_image,image,z_start,z_stop,x_start,x_stop);


%DAS_regression_top = polyval(polyfit(x_axis,meanLines.all{1},order),x_axis);
%%
figure(777);
plot(meanLines.all{8});

%%
%Check for -inf
for i = 1:length(image.all)
    if sum(isinf(meanLines.all{7})) > 2
        error('We have more than 1 inf value in the mean lateral line');
    else 
        meanLines.all{7}(1012) = meanLines.all{7}(1013);
    end
end
%%
for i = 1:length(image.all)
    regresion_coeff{i} = polyfit(x_axis,meanLines.all{i},3);
    regression_top{i} = polyval(regresion_coeff{i},x_axis);
end

for i = 1:length(image.all)
    correcting_coeff{i} = polyfit(regression_top{i},-15+linspace(0,-80,length(x_axis)),order);
    %correcting_coeff{i} = polyfit(regression_top{i},linspace(0,-80,length(x_axis)),order);
    %correcting_coeff{i} = polyfit(regression_top{i},regression_top{1},order);

end





if order == 3
    %Create calibration for "signal"
    % linear space
    x=logspace(-100/20,0,200);
    
    % dB space
    x_dB=20*log10(x);
    %x_dB_compressed=-a*x_dB.^2+b*x_dB;
    for i = 1:length(image.all)
        x_dB_compressed=correcting_coeff{i}(1)*x_dB.^3+correcting_coeff{i}(2)*x_dB.^2+correcting_coeff{i}(3)*x_dB + correcting_coeff{i}(4);
        % find the cublic spline that approximate the compressed values
        x_compressed=10.^(x_dB_compressed/20);
        gamma = fit(x.',x_compressed.','cubicspline');
        
        image_corrected.all{i} = image.all{i}.^3*correcting_coeff{i}(1) + image.all{i}.^2*correcting_coeff{i}(2) + image.all{i}*correcting_coeff{i}(3) + correcting_coeff{i}(4);
        image_corrected.all_signal{i} = reshape(gamma(abs(image.all_signal{i})./max(max(abs(image.all_signal{i})))),size(image.all_signal{i},1),size(image.all_signal{i},2));
        %image_corrected.all_signal{i} = image.all_signal{i}.^3*10.^(correcting_coeff{i}(1)/20) + image.all_signal{i}.^2*10.^(correcting_coeff{i}(2)/20) + image.all_signal{i}*10.^(correcting_coeff{i}(3)/20) + 10.^(correcting_coeff{i}(4)/20);
        
        %image_corrected.all{i}(image.all{i}<-80) = -80;
        %%
        f8888 = figure(8888);clf;
        subplot(1,2,2);
        plot(x,x,'k','LineWidth',2); hold on; grid on; axis equal tight;
        plot(x,x_compressed,'b','LineWidth',2); hold on;
        plot(x,gamma(x),'r:','LineWidth',2);
        title('Linear space');
        xlabel('Input signal');
        ylabel('Output signal');
        legend('location','nw','Uniform','p(b) mapped to linear','v(b)');
        
        %%
        f8889 = figure(8889);clf;
        subplot(1,2,1);hold all;
        plot(x_dB,x_dB,'k','LineWidth',2); hold on; grid on; axis equal tight;
        plot(x_dB,x_dB_compressed,'b','LineWidth',2); hold on;
        plot(x_dB,20*log10(gamma(x)),'r:','LineWidth',2); hold on;
        %title('Log space');
        xlabel('Input signal [dB]');
        ylabel('Output signal [dB]');
        legend('location','nw','Uniform','p(b)','20log_{10}(v(b))');
        
        f8899 = figure(8899);clf;
        subplot(1,2,1);hold all;
        plot(x_dB,x_dB,'k','LineWidth',2); hold on; grid on; axis equal tight;
        plot(x_dB,x_dB_compressed,'b','LineWidth',2); hold on;
        %title('Log space');
        xlabel('Input signal [dB]');
        ylabel('Output signal [dB]');
        legend('location','nw','Uniform','p(B)');
        xlim([-60 0]);
        
    end
   
elseif order == 2
    for i = 1:length(image.all)
        image_corrected.all{i} = image.all{i}.^2*correcting_coeff{i}(1) + image.all{i}*correcting_coeff{i}(2) + correcting_coeff{i}(3);
        %image_corrected.all{i}(image.all{i}<-80) = -80;
    end
    
elseif order == 4
    for i = 1:length(image.all)
        image_corrected.all{i} = image.all{i}.^4*correcting_coeff{i}(1) + image.all{i}.^3*correcting_coeff{i}(2) + image.all{i}.^2*correcting_coeff{i}(3) + image.all{i}*correcting_coeff{i}(4)+ correcting_coeff{i}(5);
        %image_corrected.all{i}(image.all{i}<-100) = -100;
        %image_corrected.all{i} = image_corrected.all{i}-max(image_corrected.all{i}(:));
    end
end

%[dummy,mv_idx] = ismember('MV',image.tags);
%regresion_coeff{mv_idx} = polyfit(x_axis,meanLines.all{mv_idx},3);
%regression_top{mv_idx} = polyval(polyfit(x_axis,meanLines.all{mv_idx},3),x_axis);
%correcting_coeff{mv_idx} = polyfit(regression_top{mv_idx},linspace(0,-60,length(x_axis)),3);
%correcting_coeff{mv_idx} = polyfit(regression_top{mv_idx},regression_top{1},3);


%image_corrected.all{mv_idx} = image.all{mv_idx}.^3*correcting_coeff{mv_idx}(1) + image.all{mv_idx}.^2*correcting_coeff{mv_idx}(2) + image.all{mv_idx}*correcting_coeff{mv_idx}(3) + correcting_coeff{mv_idx}(4);

% [dummy,mv_idx] = ismember('EBMV',image.tags);
% regression_top{mv_idx} = polyval(polyfit(x_axis,meanLines.all{mv_idx},3),x_axis);
% correcting_coeff{mv_idx} = polyfit(regression_top{mv_idx},regression_top{1},3);
% image_corrected.all{mv_idx} = image.all{mv_idx}.^3*correcting_coeff{mv_idx}(1) + image.all{mv_idx}.^2*correcting_coeff{mv_idx}(2) + image.all{mv_idx}*correcting_coeff{mv_idx}(3) + correcting_coeff{mv_idx}(4);


image_corrected.tags = image.tags;

%Get the mean calibrated lines;
[meanLines_full,x_axis] = getMeanLateralLines(sta_image,image_corrected,z_start,z_stop,x_start,x_stop);
[meanLines_full_calibrated,x_axis] = getMeanLateralLines(sta_image,image_corrected,z_start,z_stop,x_start,x_stop);

%%
for i = 1:length(meanLines.all)
   error_values(i) = sum(abs(meanLines_full.all{1}-meanLines_full_calibrated.all{i})); 
end

%%
figure; hold all;
for i = 1:length(regression_top)
    subplot(211);hold all
    plot(meanLines.all{i},'DisplayName',image.tags{i})
    ax(1) = gca;
    subplot(212);hold all
    plot(regression_top{i},'DisplayName',image.tags{i});
    ax(2) = gca;
end
linkaxes(ax);
legend show

%% All images 
theory=linspace(0,-80,length(x_axis));
for i = 1:length(regression_top)
    figure; 
    plot(x_axis,theory,'k--'); hold on; grid on;
    plot(x_axis,meanLines.all{i}+15,'linewidth',2);
    plot(x_axis,regression_top{i}+15);
    plot(x_axis,meanLines_full_calibrated.all{i}+15,'linewidth',2);
    xlabel('x [mm]')
    ylabel('Intensity [dB]')
    legend('Theory','Measure','Fit','Calibrated')
    title(image.tags{i});
    xlim([-20 20])
    ylim([-150 0])
end


end

