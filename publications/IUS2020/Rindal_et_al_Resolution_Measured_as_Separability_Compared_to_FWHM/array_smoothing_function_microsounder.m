
clear all;
close all;
%filename = 'P4_2v_DW134847.uff';
%filename='L7_FI_IUS2018.uff';

%FI_channel_data = uff.read_object([data_path filesep filename],'/channel_data');

%%

%FI_channel_data.probe.pitch*10^3
%FI_channel_data.probe.element_height = 7/1000;
%%
lambda_division = [1 2 4]
pitch_signs = {'\lambda','\lambda/2','\lambda/4'};
Ns = [4 8 16 32];
for i = 1:length(Ns)
    for j = 1:length(lambda_division)
        c0                      = 343;           %Speed of sound in air at 20*
        f0                      = 100*10^3;      % Center frequency [Hz]
        lambda                  = c0/f0;         % Wavelength [m]
        
        probe = uff.linear_array()
        pitch_sign = pitch_signs{j}
        probe.pitch = lambda/lambda_division(j)      %Pitch
        probe.N = Ns(i);                           %Number of elemens
        
        
        probe.element_height = probe.pitch ;
        
        %probe = FI_channel_data.probe;
        figure
        ax = subplot(211)
        probe.plot(ax)
        
        
        %%
        %lambda = FI_channel_data.lambda;
        
        
        d = probe.pitch;
        D = probe.N_elements * d;
        elemt_size = probe.element_width
        
        %Creating Kx wector
        kx_old = linspace(-3*pi/d,3*pi/d,3000);
        
        theta = linspace(-pi,pi,3000);
        kx = -(2*pi*sin(theta)/lambda);
        
        figure;hold all;
        plot(kx_old)
        plot(kx)
        
        %%
        w_m = ones(1,probe.N_elements);
        % Calculate the beam pattern
        % Input :
        %       kx  : wavenumber
        %       x   : element position
        %       w_m : element weighting
        W = zeros(length(kx),1);
        for t = 1:length(kx)
            W(t) = sum(w_m.*exp(-1i*kx(t).*probe.x'));
        end
        W = W./max(W);
        R = D;
        W_el =  sin(kx.*(elemt_size/2))./(kx/2);
        k = 2*pi/lambda;
        %W_el = sin((k.*elemt_size.*sin(theta))/2)./(k.*sin(theta)./2);
        %W_el = lambda.*sin(2*pi.*sin(theta))./(pi.*sin(theta))
        W_el = W_el./max(W_el);
        W_total = W_el.*W';
        W_total = W_total./max(W_total);
        
        figure;
        subplot(211)
        plot(db(W_el))
        
        %%
        figure;
        subplot(211);hold all;
        plot(rad2deg(theta),db(W_el))
        plot(rad2deg(theta),db(W))
        subplot(212);hold all;
        plot(rad2deg(theta),db(W_total))
        
        %%
        f = figure(123);clf;
        ax = subplot(211);hold all;
        probe.plot(ax)
        title(['Array (N=',num2str(probe.N),', width=',num2str(D*10,'%.2f'),'cm, pitch=',pitch_sign,')'],'Interpreter','tex')
        subplot(212);hold all;
        plot(rad2deg(theta),db(W),'LineWidth',1)
        plot([-90 -90],[-50 0],'r--','LineWidth',2)
        plot([90 90],[-50 0],'r--','LineWidth',2)
        ylim([-50 0])
        ylabel('Amplitude [dB]');
        xlabel('steering angle');
        title('Beampattern');
        set(gca,'FontSize',15)
        filename = ['Array_N=',num2str(probe.N),'_',strrep(pitch_sign(2:end),'/','_'),'.png']
        saveas(f,filename);
        %%
    end
end
if 0
    x_axis_mm = R.*asin(theta);
    res = R*(1.2*lambda/(D))*10^3
    res_two_way = R*(1.2*lambda/(sqrt(2)*D))*10^3
    
    [~,est_6db] = min(abs(db(W_total)+6))
    [~,est_6db_two_way] = min(abs(db(W_total.^2)+6))
    
    f = figure;clf;
    ax = subplot(311);hold all;
    %probe.plot(ax,'Probe Geometry')
    set(gca,'FontSize',14);
    subplot(312);hold all;
    plot(x_axis_mm*10^3,db(W_total),'LineWidth',2,'DisplayName','Beam Pattern');
    plot([x_axis_mm(est_6db)*10^3 x_axis_mm(est_6db)*10^3],[-60 0],'r','LineWidth',2,'DisplayName','-6 dB');
    plot([-x_axis_mm(est_6db)*10^3 -x_axis_mm(est_6db)*10^3],[-60 0],'r','LineWidth',2,'HandleVisibility','off');
    plot([-res/2 -res/2],[-60 0],'k--','LineWidth',2,'DisplayName','Res Approx Formula')
    plot([res/2 res/2],[-60 0],'k--','LineWidth',2,'HandleVisibility','off')
    set(gca,'FontSize',14);
    xlim([-2 2]);title('One Way beam pattern');xlabel('x [mm]')
    ylim([-100 0])
    legend('show','Location','se');
    subplot(313);hold all;
    plot(x_axis_mm*10^3,db(W_total.^2),'LineWidth',2,'DisplayName','Beam Pattern');
    plot([x_axis_mm(est_6db_two_way)*10^3 x_axis_mm(est_6db_two_way)*10^3],[-160 0],'r','LineWidth',2,'DisplayName','-6 dB');
    plot([-x_axis_mm(est_6db_two_way)*10^3 -x_axis_mm(est_6db_two_way)*10^3],[-160 0],'r','LineWidth',2,'HandleVisibility','off');
    plot([-res_two_way/2 -res_two_way/2],[-160 0],'k--','LineWidth',2,'DisplayName','Res Approx Formula')
    plot([res_two_way/2 res_two_way/2],[-160 0],'k--','LineWidth',2,'HandleVisibility','off')
    xlim([-2 2]);title('Two Way beam pattern');xlabel('x [mm]')
    ylim([-100 0])
    legend('show','Location','se');
    set(gca,'FontSize',14);
    set(gcf,'Position',[-1265 214 890 590]);
    saveas(f,[ustb_path,filesep,'publications',filesep,'Rindal_phd_kappe',filesep,'figures',filesep,'resolution',filesep,'beampattern_L7-4'],'eps2c')
    
    
    %%
    figure;hold all;
    plot(x_axis_mm*10^3,W_total,'DisplayName','Beam Pattern');
    plot([-res/2 -res/2],[0 1],'k--','DisplayName','Res Approx Formula')
    plot([res/2 res/2],[0 1],'k--','HandleVisibility','off')
    %plot([x_axis_mm(est_6db)*10^3 x_axis_mm(est_6db)*10^3],[-60 0],'r','LineWidth',2,'DisplayName','-6 dB');
    %plot([-x_axis_mm(est_6db)*10^3 -x_axis_mm(est_6db)*10^3],[-60 0],'r','LineWidth',2,'HandleVisibility','off');
    xlim([-1 1]);title('One Way beam pattern');xlabel('x [mm]')
end