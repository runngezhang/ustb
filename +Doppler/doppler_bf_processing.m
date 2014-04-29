clear variables;
%daq_data = '/home/inphase/Dropbox/DAQData';
%load(fullfile(daq_data,'texosetup.mat'));
path = 'F:\DAQ\packet-doppler 25.01.2012\';
daq_data = [path '2012.01.25 - packet-doppler'];
load(fullfile(path,'texosetup.mat'));
imgs = 1;

% % 
for kk=1:imgs    
    frames = (1:183) + (kk-1)*183;
    [hdr rf_ch] = load_daq_data(daq_data,ones(1,128),frames,true);  
    rf_ch = single(rf_ch);
    length(frames)
    save(fullfile(daq_data,sprintf('matlab_data_%d.mat',kk)));
end
load(fullfile(daq_data,sprintf('matlab_data_%d.mat',1)));

%%
f_demod_bmode = 7e6;
f_demod_doppler = 5e6;

iq_dec = 4;

p.fs_in = fs_rx;

p.c = 1540;         
p.dx = pitch;          
p.t0 = 2e-6;        

p.FN = 0.7;                
p.tx_angle = 0; 
p.rx_angle = 0; 

p.iq_dec = iq_dec;
p.nch = 128;

pB = p;
pD = p;

pB.fs_out = pB.fs_in/iq_dec;
pD.fs_out = pD.fs_in;

pB.f_demod = f_demod_bmode;
pD.f_demod = f_demod_doppler;

Ts_outB = 1/pB.fs_out;
Ts_outD = 1/pD.fs_out;

z_max = (pB.c/2/pB.fs_in)*size(rf_ch,1);

x = single(pB.dx*(0:0.6:(pB.nch-1)));
z_bmode = single(pB.c*pB.t0/2 + (0:pB.c*Ts_outB/2:z_max));
z_doppler = single(pD.c*pD.t0/2 + (0:pD.c*Ts_outD/2:z_max));

%%
frames_pr_frame = 83 + 100;
imgs = 1;

pB.tx_angle = tx_angles_bmode;

for kk=1:imgs        
    kk
    load(fullfile(daq_data,sprintf('matlab_data_%d.mat',kk)));
    
    iq_ch_bmode = demod_ustb(rf_ch(:,:,1:83),pB.fs_in,pB.f_demod,'narrowband',35,2e6);
    
    iq_ch_doppler = demod_ustb(rf_ch(:,:,83+(1:100)),pD.fs_in,pD.f_demod,'narrowband',35,1.5e6);
        
    [bf_dataB X_B Z_B] = planewave_beamforming2(single(iq_ch_bmode),pB,x,z_bmode,1);
    
    
    pD.tx_angle = tx_angles_doppler(1);
    pD.rx_angle = tx_angles_doppler(1);
    [iq_doppler1 X_D1 Z_D1] = planewave_beamforming2(single(iq_ch_doppler(:,:,1:2:100)),pD,x,z_doppler,1);
    
    pD.tx_angle = tx_angles_doppler(2);
    pD.rx_angle = tx_angles_doppler(2);
    [iq_doppler2 X_D2 Z_D2] = planewave_beamforming2(single(iq_ch_doppler(:,:,(2:2:100))),pD,x,z_doppler,1);
    
    iq_bmode = sum(bf_dataB,3);
    save(fullfile(daq_data,sprintf('beamformed_data_%d.mat',kk)),'iq_doppler*','iq_bmode','pD','pB','p','X_*','Z_*','tx_angles_doppler');
end

%%
env = abs(iq_bmode);
img = imagelog2(env.^2,-63,60,255);
figure(1)
clf
image(img);
colormap(gray(256))

%% Decimate doppeler data
clear variables;
daq_data = 'F:\DAQ\packet-doppler 25.01.2012\';
load(fullfile(daq_data,sprintf('beamformed_data_%d.mat',1)));
load(fullfile(daq_data,'texosetup.mat'));

env = abs(iq_bmode);
img_b = imagelog2(env.^2,-60,60,255);
figure(1)
clf
image(img_b);
colormap(gray(256))

iq_bf1 = (iq_doppler1(1:4:end-4,:,:) + iq_doppler1(2:4:end-3,:,:) + iq_doppler1(3:4:end-2,:,:) + iq_doppler1(4:4:end-1,:,:))/4;
iq_bf2 = (iq_doppler2(1:4:end-4,:,:) + iq_doppler2(2:4:end-3,:,:) + iq_doppler2(3:4:end-2,:,:) + iq_doppler2(4:4:end-1,:,:))/4;

pD.fs_out = pD.fs_out/4;

X1 = X_D1(1:4:(end-1),:);
X2 = X_D2(1:4:(end-1),:);
Z1 = Z_D1(1:4:(end-1),:);
Z2 = Z_D2(1:4:(end-1),:);

p.angles = tx_angles_doppler;
p.packetSize = 50;
p.Nstart = 1;
p.PRF = 1/(1e-9*rxD.customLineDuration)/2;
p.polyOrder = 3;
p.typeFilter = 'polyreg';
p.f_demod = pD.f_demod;
p.vNyq = p.c*p.PRF/(4*pD.f_demod);

iq_bf1 = permute(iq_bf1(:,:,p.Nstart:p.Nstart+p.packetSize-1),[3,1,2]);
iq_bf2 = permute(iq_bf2(:,:,p.Nstart:p.Nstart+p.packetSize-1),[3,1,2]);

%%
%Returns high pass filtered iq data, using polynomial regression filter (orders(0 - p-1))
[iqhp1 p.vThresh] = hp(iq_bf1,p); 
[iqhp2 p.vThresh] = hp(iq_bf2,p);

R0m1 = squeeze(mean(abs(iq_bf1).^2)); % For use with threshold trick
R11 = squeeze(mean(conj(iqhp1(2:end-1,:,:)).*iqhp1(3:end,:,:))); % R1 estimate

R0m2 = squeeze(mean(abs(iq_bf2).^2));
R12 = squeeze(mean(conj(iqhp2(1:end-1,:,:)).*iqhp2(2:end,:,:)));  

Rc1 = 0.00001*max(filter2(ones(5),R0m1)-1000,0);%+4000;
Rc2 = 0.00001*max(filter2(ones(5),R0m2)-1000,0);%+4000;

eval colormapL;
figure(2);
clf;
subplot(2,2,1);
imagesc(X1(1,:),Z1(:,1),angle(filter2(ones(3)/3,R11)+Rc1));
colormap(map1);
caxis([-pi,pi]);
colorbar; axis equal tight; 
title('Thresholding trick to exclude b-mode signal')

subplot(2,2,2);
imagesc(X1(1,:),Z1(:,1),angle(filter2(ones(3)/3,R11)));
colormap(map1);
caxis([-pi,pi]);
colorbar; 
axis equal tight; 
title('Regular phase estimate')

subplot(2,2,3);
imagesc(X2(1,:),Z2(:,1),angle(filter2(ones(3)/3,R12)+Rc2));
colormap(map1);
caxis([-pi,pi]);
colorbar; axis equal tight; 
title('Thresholding trick to exclude b-mode signal')

subplot(2,2,4);
imagesc(X2(1,:),Z2(:,1),angle(filter2(ones(3)/3,R12)));
colormap(map1);
caxis([-pi,pi]);
colorbar; 
axis equal tight; 
title('Regular phase estimate')

amphp1 = squeeze(10*log10(mean(abs(iqhp1).^2)));
amphp2 = squeeze(10*log10(mean(abs(iqhp2).^2)));
Ngray = 256;
gain = -22.5;dyn=12;

figure(3);
clf
subplot(211)
image(X1(1,:),Z1(:,1),(amphp1+gain)/dyn*Ngray);
colormap(gray(Ngray));
colorbar;
axis equal tight; 
title('IQ data after high pass filter')

subplot(212)
image(X2(1,:),Z2(:,1),(amphp2+gain)/dyn*Ngray);
colormap(gray(Ngray));
colorbar;
axis equal tight; 
title('IQ data after high pass filter')
%%
% Phase estimates for each angles, on grids 1 and 2
nrR = 6; nrB = 3;
meanFilter = ones(nrR,nrB)/(nrR*nrB);
phase1 = double(angle(filter2(meanFilter,R11)));
phase2 = double(angle(filter2(meanFilter,R12)));
phaseT1 = double(angle(filter2(meanFilter,R11)+Rc1));
phaseT2 = double(angle(filter2(meanFilter,R12)+Rc2));

X1 = double(X1);
Z1 = double(Z1);
X3 = double(X2);
Z3 = double(Z2);

aMax = max(abs(p.angles));
zMax = max(max(Z1(:)),max(Z3(:))); % Or another limit in radial direction
zMin = min(Z1(:));
xMin = 0 + zMax*tan(aMax*pi/180);
xMax = 128*304e-6 - zMax*tan(aMax*pi/180);

z = linspace(zMin,zMax,375);
x = linspace(xMin,xMax,128);

[X,Z] = meshgrid(x,z);

F = TriScatteredInterp(X1(:),Z1(:),phase1(:));
phase1 = F(X,Z);
F = TriScatteredInterp(X1(:),Z1(:),phaseT1(:));
phaseT1 = F(X,Z);
R0hp1 = double(mean(abs(iqhp1)).^2);
F = TriScatteredInterp(X1(:),Z1(:),R0hp1(:));
R0hp1 = F(X,Z);
F = TriScatteredInterp(X3(:),Z3(:),phase2(:));
phase2 = F(X,Z);
F = TriScatteredInterp(X3(:),Z3(:),phaseT2(:));
phaseT2 = F(X,Z);
R0hp2 = double(mean(abs(iqhp2)).^2);
F = TriScatteredInterp(X3(:),Z3(:),R0hp2(:));
R0hp2 = F(X,Z);

phase1(isnan(phase1)) = 0;
phase2(isnan(phase2)) = 0;
phaseT1(isnan(phaseT1)) = 0;
phaseT2(isnan(phaseT2)) = 0;
R0hp1(isnan(R0hp1)) = 0;
R0hp2(isnan(R0hp2)) = 0;

img_b_sc = interp2(X_B,Z_B,img_b,X,Z);
img_b_sc(isnan(img_b_sc)) = 0;

figure(4);
clf
subplot(2,2,1);
imagesc(X1(1,:),Z1(:,1),phaseT1);
colormap(map1);
caxis([-pi,pi]);
colorbar; 
axis equal tight; 
title('Thresholding trick to exclude b-mode signal')

subplot(2,2,2);
imagesc(X1(1,:),Z1(:,1),phase1);
colormap(map1);
caxis([-pi,pi]);
colorbar; 
axis equal tight; 
title('Regular phase estimate')

subplot(2,2,3);
imagesc(X3(1,:),Z3(:,1),phaseT2);
colormap(map1);
caxis([-pi,pi]);
colorbar; 
axis equal tight; 
title('Thresholding trick to exclude b-mode signal')

subplot(2,2,4);
imagesc(X3(1,:),Z3(:,1),phase2);
colormap(map1);
caxis([-pi,pi]);
colorbar; 
axis equal tight; 
title('Regular phase estimate')

phi = (p.angles(1)-p.angles(2) );

R0hp1(R0hp1<0) = 0;
R0hp2(R0hp2<0) = 0;

dBhp1 = 10*log10(R0hp1) - max(10*log10(R0hp1(:)));
dBhp2 = 10*log10(R0hp2) - max(10*log10(R0hp2(:)));
Pthresh = -20;
IP1 = find(dBhp1 < Pthresh);
Pthresh = -20;
IP2 = find(dBhp2 < Pthresh);

pixMat = ones(size(phase1));
pixMat(IP1) = 0;
pixMat(IP2) = 0;

figure(5);
clf; 
imagesc(pixMat);

phase1(IP1) = 0; phase2(IP1) = 0;
phase2(IP2) = 0; phase1(IP2) = 0;

phaseT1(IP1) = 0; phaseT2(IP1) = 0;
phaseT2(IP2) = 0; phaseT1(IP2) = 0;
%%
figure(6); 
clf; 
subplot(2,1,1), 
imagesc(phase1); 
colormap(map1); 
colorbar; 
axis ij; 
title('1st direction')

subplot(2,1,2), 
imagesc(phase2); 
colormap(map1); 
colorbar; 
axis ij; 
title('2nd direction'); 
set(gca,'Clim',[-pi pi]);

vT1 = phaseT1*p.PRF*p.c/(4*pi*p.f_demod);  
vT2 = phaseT2*p.PRF*p.c/(4*pi*p.f_demod);

% 
I1 = find(abs(vT1) < p.vThresh/500); %vThresh in cm/s
I2 = find(abs(vT2) < p.vThresh/500);

ixMat = ones(size(phase1));
ixMat(I1) = 0;
ixMat(I2) = 0;

figure(7); 
subplot(211)
imagesc(abs(vT1) < p.vThresh/400);
subplot(212)
imagesc(abs(vT2) < p.vThresh/400);

tixMat = pixMat+ixMat;
figure(8); 
clf; 
imagesc(tixMat);

phase1(I1) = 0;% phase2(I1) = 0;
phase2(I2) = 0;% phase1(I2) = 0;
phaseT1(I1) = 0;% phaseT2(I1) = 0;
phaseT2(I2) = 0;% phaseT1(I2) = 0;

fd1 = phaseT1*p.PRF/(2*pi);  
fd2 = phaseT2*p.PRF/(2*pi);

figure(9); 
clf; 

subplot(1,2,1), imagesc(phaseT1,[-pi pi]); 
colormap(map1); colorbar; axis ij; title('1st direction')

subplot(1,2,2), imagesc(phaseT2,[-pi pi]); 
colormap(map1); colorbar; axis ij; title('2nd direction')

% Median filtering of doppler frequency estimates
% fd1 = medfilt2(fd1,[3,3]); 
% fd2 = medfilt2(fd2,[3,3]);

%Velocity estimates, both scans:
v1 = fd1*p.c/(2*p.f_demod);
v2 = fd2*p.c/(2*p.f_demod);

% Total velocity estimate
v  = (p.c/(2*p.f_demod)) * (1/sin(phi))*sqrt( fd1.^2 + fd2.^2 - 2*fd1.*fd2*cos(phi) ); % Dunmire 2000

vx = -(p.c/(2*p.f_demod))*((fd1-fd2)./(2*sin(phi/2))); %As used in simulation scripts
vz = -(p.c/(2*p.f_demod))*((fd1+fd2)./(2*cos(phi/2))); %As used in simulation scripts (opposite sign of Kripfgans 2006 - different def of positive coordinate system and positiv velocity?)

vx = medfilt2(vx,[3,3]); 
vz = medfilt2(vz,[3,3]);

z = unique(Z(:));
x = unique(X(:));
zi = z(1:10:end);
xi = x(1:6:end);

[XI,ZI] = meshgrid(xi,zi);
vxI = interp2(X,Z,vx,XI,ZI);
vzI = interp2(X,Z,vz,XI,ZI);

LR = sqrt(vxI.^2 +vzI.^2);

figure(10); clf; pcolor(X,Z,v1); shading interp; colormap(map1); colorbar; axis ij; %set(gca,'Xlim',[0.0044 0.034])
axis equal tight; set(gca,'Clim', [-p.vNyq p.vNyq])
hold on;
quiver(XI,ZI,vxI,vzI,3,'w')
%quiver(XI,ZI,vxI./LR,vzI./LR,1,'w')
hold off;

%%
oimg = overlayImg(img_b_sc,abs(v));
oimg.setBackMap(gray(256));
oimg.setOverMap(map1);
oimg.setOverRange([-p.vNyq p.vNyq]);
oimg.setTranspRange([0 p.vNyq]);
oimg.setTranspMask(ixMat);
oimg.setAlpha(0.7);
oimg.setAxes(X(1,:),Z(:,1),X(1,:),Z(:,1));

hold(oimg.mAxHdl,'on');
quiver(oimg.mAxHdl,XI,ZI,vxI,vzI,3,'w')
hold(oimg.mAxHdl,'off');