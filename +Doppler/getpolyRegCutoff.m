function vco = getpolyRegCutoff(p)

import CoherentCompounding.Doppler.*

% General settings
    Np1 = p.packetSize;
    PRF1 = p.PRF;
    f0=p.f_demod;
    filterOrder1 = p.polyOrder;
    c = p.c;
    Nfft = 4096;

    
    %vco1 = 0.42*(filterOrder1+1)*PRF1*c/(Np1*2*f0);
   

    % Calculate polyregfilter
    FmP1 = clutterlib('poly','filtermatrix',Np1,filterOrder1);
    VP1 = clutterlib('poly','base',Np1,filterOrder1);
    H1_1 = clutterlib('poly','freqrespons',Np1,filterOrder1,Nfft,Np1,10000,1);     
    Gpoly1 = clutterlib('poly','freqrespons',Np1,filterOrder1,Nfft,Np1,10000,0); 
    GpolyLog1 = 10*log10(Gpoly1);
    vaxis1 = (freqspace(Nfft,'whole')-1)*c*PRF1/(4*f0)*100;
    vBias1 = 1540*angle(H1_1)*PRF1/(4*pi*f0)*100;    
    vNy1 = 1540*PRF1/(4*f0)*100;
 

% %% Plot results
%     figure(100);
%     hold on;
%     subplot(2,1,1);
%     set(gca,'FontSize',12,'LineWidth',1);
%     hold on;
%     plot(vaxis1,GpolyLog1,'k-','LineWidth',2);
%     xlim([-10 10]);
%     ylim([-80 5]);
%     hold off;
%     xlabel('Velocity [cm/s]');
%     ylabel('Power [dB]');
%     grid on;
%     box on;
%     title('Filter frequency responses','FontWeight','bold');
% %     legend('1','2','3','4','5','6');
% 
%     %figure(100);
%     %clf;
%     subplot(2,1,2);
%     set(gca,'FontSize',12,'LineWidth',1);
%     hold on;
%     plot(vaxis1,vBias1,'k-','LineWidth',2);
%     xlim([-10 10]);
%   
%     xlabel('Velocity [cm/s]');
%     %ylabel('Frequency bias');
%     ylabel('Velocity [cm/s]');
%     ylim([-0.7 0.7]);
%     hold off;
%     grid on;
%     box on;
%     title('Velocity bias due to clutter filter','FontWeight','bold');
% %     legend('1','2','3','4','5','6')
% %    title(['Velocity bias due to clutter filter (polyorder=' int2str(filterOrder) ', f0=' int2str(f0/1e6) 'MHz)']);
%     %legend('N60/4kHz','N12/0.8kHz','N12/4kHz','Location','SouthOutside','Orientation','horizontal');
%     
    for k = 1:1
        eval([ '[mini, ix] = min(abs(GpolyLog' num2str(k) '+3));'])
        eval(['vco(k) = abs(vaxis' num2str(k) '(ix));']);
    end
    
    