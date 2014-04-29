%% Make ROI
data_path = '/home/ingvilek/Data/Hifps Interleaved packet/';
data_path = '/home/ingvilek/Data/Patient2_2/'
%data_path = '/home/ingvilek/Data/packet-dopplerPRF4000'
data_path = '/home/ingvilek/Data/Patient2_2/'
addpath(genpath('/home/ingvilek/BeamformingCode'))

%% TO DO: Segmentation of bloodpool
% ix = 1;
% while exist(fullfile(data_path,sprintf('BloodInterpdata_%d.mat',ix)),'file') == 2
%     load( fullfile(data_path,sprintf('BloodInterpdata_%d.mat',ix)) );
%     load( fullfile(data_path,sprintf('BloodVelocitydata_%d.mat',ix)) );
%     
%     figure();
%     imagesc(v1);
%     figure();
%     imagesc(v2)
%     
%     R01 = squeeze(mean(iq_interp_1.*conj(iq_interp_1)));
%     R02 = squeeze(mean(iq_interp_2.*conj(iq_interp_2)));
%     figure(); imagesc(log10(R01))
%     figure(); imagesc(log10(R02))
%     
%     ix = ix+1;
% end

%% temporary manual segmentation                       
%rSamps = 150:320; %Ingvilds carotid
%bSamps = 50:500; % Ingvilds carotid
rSamps = 140:340; %Patient2_2
bSamps = 50:500; %Patient2_2
avgAreaSizeR = 20%20; %number of samples used to form one spectral estimate in range direction
avgAreaSizeB = 20%20; %number of samples used to form one spectral estimate in beam direction
r = 1:3:length(rSamps)-avgAreaSizeR; b = 1:5:length(bSamps)-avgAreaSizeB; % Will give overlapping submatrices 
[B R] = meshgrid(b,r); %Indeces for picking out submatrices


%% Parameters
s.nSig = 3;
s.diagonalLoading = 0;
s.setSig =0;
s.fbavg = 1;
s.clutter =0; %1--> Clutter is present
s.estSigMethod = 'energyPercentage'; % 'energyPercentage','energyProjection'
Nf = 256;

Nframes = 15;
startFrame = 34;

PRF = 4000;
spectralMatNew = zeros(Nf,length(r),length(b),Nframes);

f.typeFilter = 'polyreg';
f.polyOrder = 3;
f.packetSize = 50;
f.PRF = PRF;
f.f_demod = 5e6;
f.c = 1540;

% Spectral shift parameters
fAx = (0:1/Nf:1-(1/Nf))*PRF-PRF/2; % Original frequency axis
wLim = 200; % Frequency threshold
Plim = 100;%40 %Power threshold
fAxMat = repmat(fAx.',1,Nframes);
tAx = 1:Nframes;


%% Spectrum estimation
if 0
counter = 1;
for ix = startFrame:startFrame+Nframes-1
    load( fullfile(data_path,sprintf('BloodInterpdata_%d.mat',ix)) );

    region = iq_interp_2(1:f.packetSize,rSamps,bSamps); % Region of interest
    if s.clutter == 0
        region = hp(region,f);
    end
    
    spectralMap = {};
    matlabpool(10)

    parfor kk = 1:length(R(:))
       subMat = region(:,R(kk):R(kk)+avgAreaSizeR-1,B(kk):B(kk)+avgAreaSizeB-1); %Henter ut kolonnevis :|, |, ..., |
       %[Ptemp, lambvec] = prinCapon_orig(Nf,subMat,s);
       Pham = hamFFT(Nf,subMat); Ptemp = 10.^(Pham/10);
       
       spectralMap{kk} = Ptemp;
    end

    matlabpool close

   % figure();for kk = 1:5:length(R(:)); plot(10*log10(real(spectralMap{kk}))); pause(); end 
    spectralMat = cell2mat(spectralMap);
    spectralMatNew(:,:,:,counter) = reshape(spectralMat,[256,length(r),length(b)]);

    counter = counter+1;
end

%Reshaping
temp = permute(spectralMatNew,[1,4,2,3]);
spectrumRemse = reshape(temp,[256,Nframes,size(R,1)*size(R,2)]);


save spectrumRemse spectrumRemse
else
    load(fullfile(data_path,'spectrumRemse.mat'))
end


%% Baselineshift estimation

% Look at energy spread for different baselineshift, choose shift giving the least spread
baselineMat = zeros(size(R,1),size(R,2));
baselinePos = [100 50 0 -50 -100];

freqAxLarge = (0:1/(Nf*2):1-(1/(Nf*2)))*2-1; %Covering Twice the Nyquist limit
 for kk = 1:length(R(:))
     for bb = 1:length(baselinePos)
        shift = baselinePos(bb);
        freqAxCalc = freqAxLarge(129-shift:129+Nf-1-shift);
        freqMatCalc = repmat(freqAxCalc.',1,Nframes);

        Ptemp = circshift(squeeze(spectrumRemse(:,:,kk)),shift);
        numerator = sum(freqMatCalc.*Ptemp,1);
        denominator = sum(Ptemp,1);

        freqMean = numerator./denominator;
        freqMeanMat = repmat(freqMean,Nf,1);

        freqBW = sum( ( (freqMatCalc-freqMeanMat).^2 ).*Ptemp,1 )./denominator;
        freqBWSum(bb) = sum(freqBW,2); 

     end
     [val,baselineIX] = min(freqBWSum);
     baselineMat(kk) = baselinePos(baselineIX);

end


%%
% figure(1); 
% for kk = 1:3:length(R(:)); 
%     viewMat = squeeze(spectrumRemse(:,:,kk));
%     imagesc( 10*log10(real(circshift(viewMat,baselineMat(kk)))) - max(10*log10(real(viewMat(:) ) )) ); 
%     colormap(gray(256)); set(gca,'Clim',[-70 0]); pause(); 
% end 

%% PSV
PSVmat = zeros(size(B));
PSVCumsummat = zeros(size(B));
BloodROI = ones(size(R,1),size(R,2)); %Use std to detect whether or not we have adequate signal for spectrum display
BWthresh = Nframes*(PRF/5)^2; %Maximum tolerated BW sum in one submatrix region
for kk = 1:size(spectrumRemse,3)
    
    shift = baselineMat(kk); %Automatic baselineshift 
    
    fAx2 = (0:1/(Nf*2):1-(1/(Nf*2)))*PRF*2-PRF; %Covering twice the Nyquist range

    ftest = fAx2(129-shift:129+256-1-shift); %Shifted frequency axis
    vAx2 = ftest*1540/5e6;
    
    Ptemp = squeeze(spectrumRemse(:,:,kk));
    Ptemp2 = Ptemp;
    Pdbtemp = 10*log10(real(Ptemp))-max(10*log10(real(Ptemp(:))));
    Ptemp(Pdbtemp<-Plim) = 0;
    Ptemp(abs(fAxMat)<wLim) = 0;
    Pshift = real(circshift(Ptemp,shift));%Real power estimates for calculation of mean frequency and such
    Pshift2 = real(circshift(Ptemp2,shift)); %For display
    fAxMat2 = repmat(ftest.',1,length(tAx)); %Matrix version of shifted frequency axis

    numerator = sum(fAxMat2.*Pshift,1);
    denominator = sum(Pshift,1);

    wMean = numerator./denominator;

    wMeanMat = repmat(wMean,length(ftest),1);
   
    wBW = sum( ( (fAxMat2-wMeanMat).^2 ).*Pshift,1 )./denominator;
   
    if sum(wBW)>=BWthresh
        BloodROI(kk) = 0;
    end

    [tempmaxW,ix] = max(abs(wMean));
    maxW(kk) = tempmaxW;
    maxWix(kk) = ix;
    wBWvec(kk) = wBW(ix);
    PSV(kk) = (tempmaxW+sqrt(wBW(ix)) )*1540/5e6/2 ; %Based on mean frequency and bandwidth

    
    [a,fMeanix] = min( abs(cumsum(Pshift2) - repmat(0.5*sum(Pshift2),Nf,1 ) )); %Returns index of mean frequency for all timepoints of current submatrix
    [a,fMinix] =  min( abs(cumsum(Pshift2) - repmat(0.1*sum(Pshift2),Nf,1 ) )); %Returns index of 30% energy threshold
    [a,fMaxix] =  min( abs(cumsum(Pshift2) - repmat(0.9*sum(Pshift2),Nf,1 ) )); %Returns index of 70% energy threshold
    
    fMeanCumsum = ftest(fMeanix);
    fMinCumsum  = ftest(fMinix);
    fMaxCumsum  = ftest(fMaxix);
    
    [tempmaxF,ix] = max(abs(fMeanCumsum));
    if fMeanCumsum(ix)<0
        PSVCumsum(kk) = fMinCumsum(ix)*1540/5e6/2;
    elseif fMeanCumsum(ix)>0
        PSVCumsum(kk) = fMaxCumsum(ix)*1540/5e6/2;
    else 
        PSV = 0;
        
    end
        
    
    
    if 0 %All energy is in the clutterband beneath 200Hz
        fMeanCumsum(denominator == 0) = NaN;          
        fMinCumsum(denominator == 0) = NaN;
        fMaxCumsum(denominator == 0) = NaN;
    end
    
%     figure(2);
%     subplot(2,1,1);
%     cla; imagesc(tAx,ftest,10*log10(Pshift2)-max(10*log10(Pshift2(:)))); hold on;
%     set(gca,'Clim',[-40 0]);
%     colorbar;
%     set(gca,'YDir','Normal')
%     errorbar(tAx,wMean,sqrt(wBW),'k')
%     plot(tAx,wMean,'k','LineWidth',4)
%     colormap(gray(256))
%     hold off;
%     
%     figure(2);
%     subplot(2,1,2)
%     cla; imagesc(tAx,ftest,10*log10(Pshift2)-max(10*log10(Pshift2(:)))); hold on;
%     set(gca,'Clim',[-40 0]);
%     colorbar;
%     set(gca,'YDir','Normal')
%     plot(tAx,fMeanCumsum,'k','LineWidth',4);
%     plot(tAx,fMinCumsum,'k:')
%     plot(tAx,fMaxCumsum,'k:')
%     colormap(gray(256))
%     hold off;
end

PSVCumsummat(:) = PSVCumsum; %Peak systolic velocity based on cumulative Power sum

PSVmat(:) = PSV;

maxWmat = zeros(size(PSVmat));
maxWmat(:)= maxW;

N = hist(maxWix,1:Nframes);
[maxHist,histix] = max(N);
systolicFrame = startFrame+histix-1;

%% Angle estimation and final figures

load( fullfile(data_path,sprintf('BloodVelocitydata_%d.mat',systolicFrame)) );
v2region = v2(rSamps,bSamps);
vxregion = vx(rSamps,bSamps);
vzregion = vz(rSamps,bSamps);
angles = zeros(size(v2region));
phi = 20*pi/180;

 
 for k = 1:length(v2region(:))                                                 
            VZ = vzregion(k); VX = vxregion(k);                                         

            if VZ<0 && VX>0                                                 
                angles(k) = pi/2 -atan(abs(VZ)/abs(VX)) +abs(phi/2);
            elseif VZ<0 && VX<0
                angles(k) = pi/2 -atan(abs(VZ)/abs(VX)) -abs(phi/2);
            elseif VZ >0 && VX <0
                angles(k) = pi/2 +atan(abs(VZ)/abs(VX)) -abs(phi/2);
            elseif VZ >0 && VX >0 
                angles(k) = pi/2 +atan(abs(VZ)/abs(VX)) +abs(phi/2);
            else
                angles(k) = NaN;
            end
 end
 
angleMat = zeros(size(R));
 
for kk = 1:length(R(:))
     subMat = angles(R(kk):R(kk)+avgAreaSizeR-1,B(kk):B(kk)+avgAreaSizeB-1); %Henter ut kolonnevis :|, |, ..., |  
     angleMat(kk) = mean(subMat(:));
end
angleMat1 = angleMat-phi; %From direction1

PSVCumsummat(BloodROI==0) = NaN;

figure(); imagesc(angleMat*180/pi); title('Angle with second txrxdirection')
figure(); imagesc(maxWmat./cos(angleMat)*1540/5e6/2); title('Mean Velocity from spectrum angle corrected'); set(gca,'Clim',[-1 1])
figure(); imagesc(PSVmat./cos(angleMat)); title('Mean velocity +SD from spectrum angle corrected'); set(gca,'Clim',[-1 1])
figure(); imagesc(PSVCumsummat./cos(angleMat)); title('PSV from cumsum'); set(gca,'Clim',[-1 1])





    