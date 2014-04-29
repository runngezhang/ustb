clear variables;
%close all;
%data_path = 'C:\Documents and Settings\ingvilek\Desktop\Jobb\Ultrasonix\';
%data_path = 'F:\DAQ\2011.12.04 - COMBO\';
data_path = 'F:\DAQ\2011.12.05 - COMBO INGVILD\';
angles = [-10 10];
type_data = 'iq';% 'rf', 'iq';

if strcmp(type_data,'iq')
    load([data_path 'iq_bf_data_' num2str(angles(1))]);
    load([data_path 'iq_bf_data_' num2str(angles(2))]);
elseif strcmp(type_data,'rf')
    load([data_path 'rf_bf_data_' num2str(angles(1))]);
    load([data_path 'rf_bf_data_' num2str(angles(2))]);
else
    disp('Data type not supported. Use rf or iq');
end
eval colormapL;

p.angles = angles;
%p.packetSize = 50;
p.Nstart = 1;
p.PRF = 4000;
p.polyOrder =3;
p.typeFilter = 'polyreg';

p.vNyq = p.c*p.PRF/(4*p.f_demod);

%% R1 estimates and threshold info
iq_bf1 = permute(iq_bf1(:,:,p.Nstart:p.Nstart+p.packetSize-1),[3,1,2]);
iq_bf2 = permute(iq_bf3(:,:,p.Nstart:p.Nstart+p.packetSize-1),[3,1,2]);

[iqhp1 p.vThresh] = hp(iq_bf1,p); %Returns high pass filtered iq data, using polynomial regression filter (orders(0 - p-1))
[iqhp2 p.vThresh] = hp(iq_bf2,p);

disp(['-3dB velocity after clutter filter is ' num2str(p.vThresh) ' cm/s' ])

R0m1=squeeze(mean(abs(iq_bf1).^2)); % For use with threshold trick
R11=squeeze(mean(conj(iqhp1(2:end-1,:,:)).*iqhp1(3:end,:,:))); % R1 estimate

R0m2=squeeze(mean(abs(iq_bf2).^2));
R12=squeeze(mean(conj(iqhp2(1:end-1,:,:)).*iqhp2(2:end,:,:)));  

%Rc1=0.00001*max(filter2(ones(5),R0m1)-1000,0)+4000; %Ask Hans %original
%Rc2=0.00001*max(filter2(ones(5),R0m2)-1000,0)+4000; %Ask Hans %original
Rc1 = 0.00001*max(filter2(ones(5),R0m1)-1000,0);%+4000;
Rc2 = 0.00001*max(filter2(ones(5),R0m2)-1000,0);%+4000;


figure(1);
clf;
subplot(2,1,1);
imagesc(X1(:),Z1(:),angle(filter2(ones(3),R12)+Rc2));colormap(map1);caxis([-pi,pi]);colorbar; axis equal tight; title('Thresholding trick to exclude b-mode signal')
subplot(2,1,2);
imagesc(X1(:),Z1(:),angle(filter2(ones(3),R12)));colormap(map1);caxis([-pi,pi]);colorbar; axis equal tight; title('Regular phase estimate')

figure(2);
clf
amphp=squeeze(10*log10(mean(abs(iqhp2).^2)));
Ngray=256;
gain=-10;dyn=12;
image(X1(:),Z1(:),(amphp+gain)/dyn*Ngray);colormap(gray(Ngray));colorbar;axis equal tight; title('IQ data after high pass filter')

% Phase estimates for each angles, on grids 1 and 2
nrR = 6; nrB = 3;
meanFilter = ones(nrR,nrB)/(nrR*nrB);
phase1 = double(angle(filter2(meanFilter,R11)));
phase2 = double(angle(filter2(meanFilter,R12)));
phaseT1 = double(angle(filter2(ones(3),R11)+Rc1));
phaseT2 = double(angle(filter2(ones(3),R12)+Rc2));

X1 = double(X1);
Z1 = double(Z1);
X3 = double(X3);
Z3 = double(Z3);


%% Interpolate phase estimates to common grid
aMax = max(abs(p.angles));
zMax = max(max(Z1(:)),max(Z3(:))); % Or another limit in radial direction
zMin = 0;
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

figure(3);
clf
subplot(2,1,1);
imagesc(X1(:),Z1(:),phaseT2);colormap(map1);caxis([-pi,pi]);colorbar; axis equal tight; title('Thresholding trick to exclude b-mode signal')
subplot(2,1,2);
imagesc(X1(:),Z1(:),phase2);colormap(map1);caxis([-pi,pi]);colorbar; axis equal tight; title('Regular phase estimate')


save(fullfile('phase_data.mat'),'phase1','phase2','p','phaseT1','phaseT2','X','Z','R0hp1','R0hp2');

%% Vector doppler
clear variables;
eval colormapL;
load phase_data.mat
phi = (p.angles(1)-p.angles(2) )*pi/180;

R0hp1(R0hp1<0) = 0;
R0hp2(R0hp2<0) = 0;

dBhp1 = 10*log10(R0hp1) - max(10*log10(R0hp1(:)));
dBhp2 = 10*log10(R0hp2) - max(10*log10(R0hp2(:)));
Pthresh = -5;
IP1 = find(dBhp1 < Pthresh);
IP2 = find(dBhp2 < Pthresh);

pixMat = ones(size(phase1));
pixMat(IP1) = 0;
pixMat(IP2) = 0;

figure(4);clf; imagesc(pixMat);

phase1(IP1) = 0; phase2(IP1) = 0;
phase2(IP2) = 0; phase1(IP2) = 0;

figure(5); clf; subplot(2,1,1), imagesc(phase1); colormap(map1); colorbar; axis ij; title('1st direction')
subplot(2,1,2), imagesc(phase2); colormap(map1); colorbar; axis ij; title('2nd direction'); set(gca,'Clim',[-pi pi]);

vT1 = phaseT1*p.PRF*p.c/(4*pi*p.f_demod);  
vT2 = phaseT2*p.PRF*p.c/(4*pi*p.f_demod);

% v1 = phase1*p.PRF*p.c/(4*pi*p.f_demod);  
% v2 = phase2*p.PRF*p.c/(4*pi*p.f_demod);

% 
I1 = find(abs(vT1) < p.vThresh/100); %vThresh in cm/s
I2 = find(abs(vT2) < p.vThresh/100);

% I1 = find(abs(v1) < p.vThresh);
% I2 = find(abs(v2) < p.vThresh);


ixMat = ones(size(phase1));
ixMat(I1) = 0;
ixMat(I2) = 0;

figure(6); clf; imagesc(ixMat);

tixMat = pixMat+ixMat;
figure(7); clf; imagesc(tixMat);

phase1(I1) = 0; phase2(I1) = 0;
phase2(I2) = 0; phase1(I2) = 0;
fd1 = phase1*p.PRF/(2*pi);  
fd2 = phase2*p.PRF/(2*pi);

figure(8); clf; subplot(2,1,1), imagesc(phase1); colormap(map1); colorbar; axis ij; title('1st direction')
subplot(2,1,2), imagesc(phase2); colormap(map1); colorbar; axis ij; title('2nd direction')

% Median filtering of doppler frequency estimates
fd1 = medfilt2(fd1,[3,3]); 
fd2 = medfilt2(fd2,[3,3]);

%Velocity estimates, both scans:
v1 = fd1*p.c/(2*p.f_demod);
v2 = fd2*p.c/(2*p.f_demod);

% Total velocity estimate
v  = (p.c/(2*p.f_demod)) * (1/sin(phi))*sqrt( fd1.^2 + fd2.^2 - 2*fd1.*fd2*cos(phi) ); % Dunmire 2000

% Velocity vector estimates
% vx = -(p.c/(2*p.f_demod))*((fd1-fd2)./(2*sin(phi/2))); %Kripfgans 2006 (Valid only for parallel tx/rx...)
% vz = +(p.c/(2*p.f_demod))*((fd1+fd2)./(2*cos(phi/2))); %Kripfgans 2006 (Valid only for parallel tx/rx...)

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

figure(9); clf; pcolor(X,Z,v1); shading interp; colormap(map1); colorbar; axis ij; %set(gca,'Xlim',[0.0044 0.034])
axis equal tight; set(gca,'Clim', [-p.vNyq p.vNyq])
hold on;
quiver(XI,ZI,vxI,vzI,3,'w')
%quiver(XI,ZI,vxI./LR,vzI./LR,1,'w')
hold off;


