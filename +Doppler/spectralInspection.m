 %% GENERATE SPECTRAL ESTIMATES
  clear all; close all;
  dataPath = 'F:\In vivo bifurcation\packet-dopplerPRF4000\'; %In vivo recording of Ingvilds carotid bifurcation
  load([dataPath 'beamformed_data_1']);
  env = abs(iq_bmode); clear iq_bmode; clear Z_B Z_D1 Z_D2 X_D1 X_D2 X_B;
  img1 = 30;
  imgs = 45;
   
  figure(300); imagesc(env);
  [xp,zp] = ginput(1);
         
  for kk = img1:imgs  
      clearvars -except kk xp zp dataPath img1 imgs;
      addR = 9;
      addB = 9;
      
      Nf = 256;
      PRF = 4000; 
      shift = 70;
      f = ((1/Nf)*(0:Nf-1)*PRF -PRF/2);
      df = f(2)-f(1);
      f = (f-df*shift)*(-1);
      v = f*1540/(2*5e6);
     
      ranges = round(zp)-addR:1:round(zp)+addR;
      beams =  round(xp)-addB:1:round(xp)+addB;
      
      load([dataPath sprintf('beamformed_data_%i',kk)]);
      clear iq_bmode; clear Z_B Z_D1 Z_D2 X_D1 X_D2 X_B;
            
      p.angles = tx_angles_doppler;
      p.packetSize = 50;
      p.PRF = PRF;% 1/(1e-9*rxD.customLineDuration)/2; %Not entirely true. 
      p.polyOrder = 3;
      p.typeFilter = 'hSTLow1';
      p.f_demod = pD.f_demod;
      p.vNyq = p.c*p.PRF/(4*pD.f_demod); 
      
      s.diagonalLoading = 0;
      s.fbavg = 1;
      s.nSig = 10;
      s.setSig = 0;
      
      iq_bf1 = permute((iq_doppler1(1:4:end-4,:,:) + iq_doppler1(2:4:end-3,:,:) + iq_doppler1(3:4:end-2,:,:) + iq_doppler1(4:4:end-1,:,:))/4 ,[3,1,2]);
      %iq_bf2 = (iq_doppler2(1:4:end-4,:,:) + iq_doppler2(2:4:end-3,:,:) + iq_doppler2(3:4:end-2,:,:) + iq_doppler2(4:4:end-1,:,:))/4;
      clear iq_doppler1;
      %clear iq_doppler2;
    
      % High pass clutter filter and compute the correlation functions
      [iqhp1 p.vThresh] = hp(iq_bf1,p); 
    %  [iqhp2 p.vThresh] = hp(iq_bf2,p);
      
      ROIham = iqhp1(:,ranges,beams);
      pwHam = hamFFT(Nf,ROIham);
      clear iqhp1;
      
      quality = 10;
      fac = floor(p.packetSize/quality);
      cutPackS = fac*quality;
      ROIcap = iq_bf1(:,ranges,beams);
      ROIcap  = ROIcap (1:cutPackS,:,:);
      [K,R,B,frames] = size(ROIcap );
      ROIcap  = permute(ROIcap ,[1,4, 2, 3]);
      ROIcap  = reshape(ROIcap ,[K/fac, frames*fac, R, B]);
      ROIcap = permute(ROIcap ,[1,3, 4, 2]); % Reshaped version

      [pwCap,lambda] = prinCapon(Nf,ROIcap,s);
      
      save([dataPath sprintf('pwCap_ensemble10_%i',kk)],'pwCap');
      save([dataPath sprintf('pwHam_ensemble45_%i',kk)],'pwHam');
  end
  
  
  
  %% VIEW SPECTRAL ESTIMATES
  
  pwCapFull = zeros(Nf,(imgs-img1+1)*fac);
  pwHamFull = zeros(Nf,(imgs-img1+1)); 
  ix = 1;
  for kk = img1:imgs
      load([dataPath sprintf('pwCap_ensemble10_%i',kk)]);
      load([dataPath sprintf('pwHam_ensemble45_%i',kk)]);

      pwCapFull(:,(ix-1)*fac +1:(ix-1)*fac+fac) = pwCap;
      pwHamFull(:,ix) = pwHam;
      ix = ix+1;
  end
  figure(); imagesc(10*log10(abs(pwCapFull)/max(abs(pwCapFull(:)))))
  set(gca,'Clim',[-60 0]);
  colormap(gray(256))
  
  figure(); imagesc(pwHamFull -max(pwHamFull(:)))
  colormap(gray(256))
    set(gca,'Clim',[-20 0]);

  
  