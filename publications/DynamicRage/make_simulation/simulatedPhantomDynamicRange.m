function [sca_tot,amp_tot] = simulatedPhantomDynamicRange(sca_per_mm2)


%% Create lateral gradient (lg)
x_min_lg = -15/1000;
x_max_lg = 15/1000;
z_min_lg = 47.5/1000;
z_max_lg = 52.5/1000;
Intensity_lg = 0;
dB_mm_lg = 2;

[sca_lg,amp_lg] = simulatedPhantomGradientBlock(sca_per_mm2,x_min_lg,x_max_lg,z_min_lg,z_max_lg,Intensity_lg,dB_mm_lg)

x_min_lg_left = -25/1000;
x_max_lg_left = -15/1000;
z_min_lg_left = 47.5/1000;
z_max_lg_left = 52.5/1000;
Intensity_lg_left = 0;
dB_mm_lg_left = 0;

[sca_lg_left,amp_lg_left] = simulatedPhantomGradientBlock(sca_per_mm2,x_min_lg_left,x_max_lg_left,z_min_lg_left,z_max_lg_left,Intensity_lg_left,dB_mm_lg_left)

x_min_lg_right = 15/1000;
x_max_lg_right = 25/1000;
z_min_lg_right = 47.5/1000;
z_max_lg_right = 52.5/1000;
Intensity_lg_right = -abs(x_min_lg-x_max_lg)*1000*dB_mm_lg;
dB_mm_lg_right = 0;
[sca_lg_right,amp_lg_right] = simulatedPhantomGradientBlock(sca_per_mm2,x_min_lg_right,x_max_lg_right,z_min_lg_right,z_max_lg_right,Intensity_lg_right,dB_mm_lg_right)


%% 

%% Create axial gradient (ag)
x_min_ag = 10/1000;
x_max_ag = 15/1000;
z_min_ag = 12.5/1000;
z_max_ag = 42.5/1000;
Intensity_ag = 0;
dB_mm_ag = 2;

[sca_ag,amp_ag] = simulatedPhantomGradientBlock(sca_per_mm2,x_min_ag,x_max_ag,z_min_ag,z_max_ag,Intensity_ag,[0 dB_mm_ag])

x_min_ag_top = 10/1000;
x_max_ag_top = 15/1000;
z_min_ag_top = 2.5/1000;
z_max_ag_top = 12.5/1000;
Intensity_ag_top = 0;
dB_mm_ag_top = 0;

[sca_ag_top,amp_ag_top] = simulatedPhantomGradientBlock(sca_per_mm2,x_min_ag_top,x_max_ag_top,z_min_ag_top,z_max_ag_top,Intensity_ag_top,[0 dB_mm_ag_top])

x_min_ag_bottom = 10/1000;
x_max_ag_bottom = 15/1000;
z_min_ag_bottom = 42.5/1000;
z_max_ag_bottom = 45.5/1000;
Intensity_ag_bottom= -abs(z_min_ag-z_max_ag)*1000*dB_mm_ag;
dB_mm_ag_bottom = 0;

[sca_ag_bottom,amp_ag_bottom] = simulatedPhantomGradientBlock(sca_per_mm2,x_min_ag_bottom,x_max_ag_bottom,z_min_ag_bottom,z_max_ag_bottom,Intensity_ag_bottom,[0 dB_mm_ag_bottom])


%% Create hypoechoic cyst (hc)
x_min_hc = -15/1000;
x_max_hc = 0/1000;
z_min_hc = 17.5/1000;
z_max_hc = 32.5/1000;
Intensity_hc = -10;
dB_mm_hc = 0;

[sca_hc_temp,amp_hc_temp] = simulatedPhantomGradientBlock(sca_per_mm2,x_min_hc,x_max_hc,z_min_hc,z_max_hc,Intensity_hc,dB_mm_hc)

% Define geometry of nonechoic cyst
r=2.5e-3;           % Radius of cyst [m]
xc=-7.5e-3;           % Position of cyst in x [m]
zc=27.5e-3;           % Position of cyst in z [m]
%Find the indexes inside cyst
inside_nonechoic = (((sca_hc_temp(:,1)-xc).^2 + (sca_hc_temp(:,3)-zc).^2) < r^2);

%%
clear sca_hc;
sca_hc(:,1) = sca_hc_temp(inside_nonechoic==0,1);
sca_hc(:,2) = sca_hc_temp(inside_nonechoic==0,2);
sca_hc(:,3) = sca_hc_temp(inside_nonechoic==0,3);
amp_hc = amp_hc_temp(inside_nonechoic==0);

%% Create amplitude box
x_min_left_l = -15/1000;
x_max_left_l = -12.5/1000;
z_min_left_l = 35/1000;
z_max_left_l = 40/1000;
Intensity_left_l = 0;

[sca_left_l,amp_left_l] = simulatedPhantomGradientBlock(sca_per_mm2,x_min_left_l,x_max_left_l,z_min_left_l,z_max_left_l,Intensity_left_l,0)
x_min_left_r = -12.5/1000;
x_max_left_r = -10/1000;
Intensity_left_r = -10;
[sca_left_r,amp_left_r] = simulatedPhantomGradientBlock(sca_per_mm2,x_min_left_r,x_max_left_r,z_min_left_l,z_max_left_l,Intensity_left_r,0)

x_min_right_l = -5/1000;
x_max_right_l = -2.5/1000;
z_min_right_l = 35/1000;
z_max_right_l = 40/1000;
Intensity_right_l = 0;

[sca_right_l,amp_right_l] = simulatedPhantomGradientBlock(sca_per_mm2,x_min_right_l,x_max_right_l,z_min_right_l,z_max_right_l,Intensity_right_l,0)
x_min_right_r = -2.5/1000;
x_max_right_r = 0/1000;
Intensity_right_r = -35;
[sca_right_r,amp_right_r] = simulatedPhantomGradientBlock(sca_per_mm2,x_min_right_r,x_max_right_r,z_min_right_l,z_max_right_l,Intensity_right_r,0)



%%
sca_points(1,:) =[-7.5e-3,  0, 15e-3];    % point scatterer position [m]
sca_points(2,:) =[-7.5e-3,  0, 20e-3];    %Pont in speckle
sca_points(3,:) =[-7.5e-3,  0, 45e-3];    % point scatterer position [m]

amp_points = max([amp_lg(:); amp_hc(:); amp_ag(:)]).*ones(length(sca_points),1);
amp_points(2) = amp_points(3)*1.2;
%%

sca_tot = [sca_lg; sca_lg_left; sca_lg_right; sca_ag_top; sca_hc; sca_ag; sca_ag_bottom; sca_left_l; sca_left_r; sca_right_l; sca_right_r; sca_points];
amp_tot = [amp_lg; amp_lg_left; amp_lg_right; amp_ag_top; amp_hc; amp_ag; amp_ag_bottom; amp_left_l; amp_left_r; amp_right_l;  amp_right_r; amp_points];

figure;
plot3(sca_tot(:,1)*1e3,sca_tot(:,3)*1e3,20*log10(amp_tot),'b.'); axis equal;
zlim([-60 5]);
xlabel('X [mm]');ylabel('Z [mm]');zlabel('Amplitude [dB]');