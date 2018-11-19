function [diff] = dynamic_range_test(channel_data,b_data,title_txt)

if strcmp(channel_data.name,'Simulated dynamic range phantom. Created with Field II. See the reference for details')
    z_start = 47.5;
    z_stop = 52.5;
    x_start = -12.5;
    x_stop = 12.5;
    gradient = -2;% db/mm
    do_axial = 2; %If we want to estimate the axial as well, run this twice.
    sub_fig_setup = [1 3];
 %   theory_lateral = -40*(b_data.scan.x_axis(mask_lateral)+10e-3)/20e-3;
elseif strcmp(channel_data.name,'New Simulated dynamic range phantom. Created with Field II. See the reference for details')
    z_start = 40;
    z_stop = 48.5;
    x_start = -12.05;
    x_stop = 12.05;
    gradient = -1.66;% db/mm
    do_axial = 2; %If we want to estimate the axial as well, run this twice.
    sub_fig_setup = [1 3];
elseif strcmp(channel_data.name,'v4 New Simulated dynamic range phantom. Created with Field II. See the reference for details')
    z_start = 40;
    z_stop = 48.5;
    x_start = -15;
    x_stop = 15;
    gradient = -2;% db/mm
    do_axial = 2; %If we want to estimate the axial as well, run this twice.
    sub_fig_setup = [1 3];
elseif strcmp(channel_data.name,'Experimental dynamic range phantom. Recorded on a Verasonics Vantage 256 with a L11 probe. See the reference for details')
    z_start = 40;
    z_stop = 48.5;
    x_start = -15;
    x_stop = 15;
    gradient = -1.66;% db/mm
    do_axial = 1;
    sub_fig_setup = [1 2];
else
    error('The dynamic range test is only defined for the simulated dynamic range phantom, and the experimental one. See the example.');
end

figure(); hold all;
for i = 1:do_axial
    % Mask out the top gradient part of the image
    mask=reshape(b_data.scan.z>z_start*10^-3&b_data.scan.z<z_stop*10^-3,[length(b_data.scan.z_axis) length(b_data.scan.x_axis)])...
        &reshape(b_data.scan.x>x_start*10^-3&b_data.scan.x<x_stop*10^-3,[length(b_data.scan.z_axis) length(b_data.scan.x_axis)]);
    
    % Find the x min and max index in the image from the gradient mask, top
    z_min = rem(min(find(mask==1)),length(b_data.scan.z_axis));
    z_max = rem(max(find(mask==1)),length(b_data.scan.z_axis));
    
    temp = find(mask(z_min,:)==1);
    x_min = temp(1);
    x_max = temp(end);
    
    img = b_data.get_image();
    mean_line = mean(img(z_min:z_max,x_min:x_max),i);
    mean_line = mean_line-max(mean_line(:)); %Normalize to have the gradient start at 0
    
    if i == 1 %Then we are estimating the lateral
        %regresion_coeff = polyfit(b_data.scan.x_axis(x_min:x_max)*10^3,lateral_line',1);
        %regression = polyval(regresion_coeff,b_data.scan.x_axis(x_min:x_max)*10^3);
        theory_line = (gradient*(b_data.scan.x_axis(x_min:x_max)-x_start*10^-3))*10^3;

        
        sub_fig_index_1 = [1];
        sub_fig_index_2 = [2];
    else % we are estimating the axial
        %regresion_coeff = polyfit(b_data.scan.z_axis(z_min:z_max)*10^3,mean_line,1);
        %regression = polyval(regresion_coeff,b_data.scan.z_axis(z_min:z_max)*10^3);
        theory_line = (gradient*(b_data.scan.z_axis(z_min:z_max)-z_start*10^-3))*10^3;
        sub_fig_index_1 = [1];
        sub_fig_index_2 = [2];
    end
    
    regresion_coeff = polyfit(theory_line,mean_line(:),1);
    regression = polyval(regresion_coeff,theory_line);
    
    if i == 1
        subplot(sub_fig_setup(1),sub_fig_setup(2),sub_fig_index_1); hold all;
        imagesc(b_data.scan.x_axis*1e3,b_data.scan.z_axis*1e3,b_data.get_image());
        colormap gray; axis image; caxis([-60 0]);
        xlabel('x [mm]');ylabel('z [mm]');
        set(gca,'FontSize',15);if(nargin == 3); title(title_txt); end
        set(gca, 'YDir','reverse');
        c = 'r';
    else
        c = 'g'; 
    end
    subplot(sub_fig_setup(1),sub_fig_setup(2),sub_fig_index_1); hold all;
    rectangle('Position',[x_start z_start abs(x_start-x_stop) abs(z_start-z_stop)],...
         'LineWidth',2,'LineStyle','--','EdgeColor',c)

    subplot(sub_fig_setup(1),sub_fig_setup(2),sub_fig_index_2+i-1);
    hold all;
    plot(theory_line,theory_line,'LineWidth',2,'DisplayName','Theoretical');
    plot(theory_line,mean_line,'LineWidth',2,'DisplayName','Mean lateral line');
    plot(theory_line,regression,'LineWidth',2,'DisplayName','Estimated slope');
    set(gca, 'XDir','reverse');xlabel('Input [dB]');ylabel('Output [dB]');
    title(sprintf('Theory: -1, estimated: -%.4f',regresion_coeff(1)));
    legend show;set(gca,'FontSize',15); ylim([-80 0]);
    
    v= max(abs(theory_line));
      rectangle('Position',[-v -80 v 80],...
         'LineWidth',4,'LineStyle','-','EdgeColor',c)
    xlim([-v 0]);
    %axis([-v 0 -v 0])
    
    diff(i) = abs(1-regresion_coeff(1))*100;
    
    if(do_axial == 2)
        if strcmp(channel_data.name,'Simulated dynamic range phantom. Created with Field II. See the reference for details')
            z_start = 12.5;
            z_stop = 40;
            x_start = 10.5;
            x_stop = 14.5;
        elseif strcmp(channel_data.name,'New Simulated dynamic range phantom. Created with Field II. See the reference for details')
            z_start = 10;
            z_stop = 40;
            x_start = 10.5;
            x_stop = 19;
        elseif strcmp(channel_data.name,'v4 New Simulated dynamic range phantom. Created with Field II. See the reference for details')
            z_start = 13;
            z_stop = 37;
            x_start = 15;
            x_stop = 18.5;
        end
    end
    
%     annotation('textbox', [0 0.5 1 0.1], ...
%     'String', 'hello, title', ...
%     'EdgeColor', 'none', ...
%     'HorizontalAlignment', 'center',...
%     'TextSize',20)
    
end

end

