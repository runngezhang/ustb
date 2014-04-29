  %% ..check after running process_pw_doppler (high pass filtered iq data should be in the form [ranges, beams, packet])
  %iq_bf1 = permute(iq_bf1,[2,3,1]);
  %iq_bf2 = permute(iq_bf2,[2,3,1]);
  %iq_bfx = iq_bf2;
  %iq_bfx = permute(iqhp1,[2,3,1]); %IQ data ARE NOT interpolated
  iq_bfx = permute(iqhp_interp_1,[2,3,1]); %IQ data ARE interpolated
  
  decimate = 1; %Decimate using this factor.
 
  iq_bfx = iq_bfx(:,:,1:decimate:end);
  
  while true %ctrl C to stop...:-)
      figure(200);
      clf;
      subplot(1,2,1);
      imagesc(X_B(:),Z_B(:),20*log10( env/max(env(:)) ) ); colormap(gray(256)); 
      set(gca,'Clim',[-60 0]);

      [xin,zin] = ginput(1);
      %xp = 20.77/1000; zp = 19.984/1000;
      xp = length(find( x < xin ));
      zp = length(find( z < zin ));
      
      
      h = impoint(gca,xin,zin);
      fcn = makeConstrainToRectFcn('impoint',get(gca,'XLim'),get(gca,'YLim'));
      setPositionConstraintFcn(h,fcn);

      
      %xp = 108.7674; zp = 246.5088;
      subplot(3,2,2)
      plot(squeeze(real(iq_bfx(round(zp),round(xp),:))),'b','DisplayName','real');
      hold on;
      plot(squeeze(imag(iq_bfx(round(zp),round(xp),:))),'r','DisplayName','imag')
      legend('show')
      xlabel('Slow time [packet]'); ylabel('Amplitude')
      hold off;

      Nfft = 128;
      PRF = 4000; %PRF = p.PRF;
      PRF = PRF/decimate;
      f = linspace(-PRF/2,PRF/2,Nfft);
      v = f*1540/(2*5e6);
      p.packetSize = size(iq_bfx,3);
      
      s.diagonalLoading = 0;
      s.fbavg = 1;
      s.nSig = 10;
      s.setSig = 10;
      addR = 9; %9x9 gives 1x1mm ROI with axes x,z
      addB = 9;
      
      fftSig = fftshift( fft(squeeze(iq_bfx(round(zp),round(xp),:)),Nfft));

      subplot(3,2,4)
      plot(f,abs(fftSig));

      subplot(3,2,6);
      plot(f,cumsum(abs(fftSig)/sum(abs(fftSig))),'DisplayName','Cumulative absolute value')
      set(gca,'Ylim',[0 1]);  xlabel('Frequency [Hz]'); ylabel('abs(fftSig)');
      hold on;

      [a,fMeanix] = min( abs(cumsum(abs(fftSig)) - 0.5*sum(abs(fftSig)) ));
      [a,fMinix] = min( abs(cumsum(abs(fftSig)) - 0.3*sum(abs(fftSig)) ));
      [a,fMaxix] = min( abs(cumsum(abs(fftSig)) - 0.7*sum(abs(fftSig)) ));

      plot([f(fMeanix),f(fMeanix)],[0 1],':r','DisplayName',sprintf('vMean = %f',v(fMeanix)));
      plot([f(fMinix),f(fMinix)],[0 1],':k','DisplayName',sprintf('vMin = %f',v(fMinix)));
      plot([f(fMaxix),f(fMaxix)],[0 1],':k','DisplayName',sprintf('vMax = %f',v(fMaxix)));

      legend('show')
      hold off;
      subplot(1,2,1); title('Press any key to continue')
    
      ranges = round(zp)-addR:1:round(zp)+addR;
      beams =  round(xp)-addB:1:round(xp)+addB;
      
%       ROIham = iqhp1(:,ranges,beams);
      ROIham = permute(iq_bfx(ranges,beams,:),[3,1,2]);
      pwHam = hamFFT(Nfft,ROIham);
      
      quality = 13;
      fac = floor(p.packetSize/quality);
      cutPackS = fac*quality;
%       ROIcap = iqhp1(:,ranges,beams);
      ROIcap = permute(iq_bfx(ranges,beams,:),[3,1,2]);
      ROIcap  = ROIcap (1:cutPackS,:,:);
      [K,R,B,frames] = size(ROIcap );
      ROIcap  = permute(ROIcap ,[1,4, 2, 3]);
      ROIcap  = reshape(ROIcap ,[K/fac, frames*fac, R, B]);
      ROIcap = permute(ROIcap ,[1,3, 4, 2]); % Reshaped version
      
      [pwCap,lambda] = prinCapon(Nfft,ROIcap,s);
      
      figure(4); clf; plot(v,pwHam - max(pwHam(:)),'k','LineWidth',2,'DisplayName', 'hamFFT on hp iq data');
      hold on;
      for ix = 1:size(pwCap,2);
          plot(v,10*log10(abs(pwCap(:,ix))/max(abs(pwCap(:,ix)))),'LineStyle',':','Color',[1/ix, ix/10, 1/ix],'DisplayName',sprintf('pCapon on hp iq data %i',ix));
      end
      
      axis tight;
      title(sprintf('Temporal resolution of pCapon is %i times of pHam. The number of signal vectors is %i.',fac,lambda))
      hold off;
      legend show;
     
      
      pause();
  end
