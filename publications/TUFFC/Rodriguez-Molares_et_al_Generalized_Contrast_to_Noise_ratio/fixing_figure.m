
fig = openfig('Final_figure_NLM.fig');

url='https://nyhirse.medisin.ntnu.no/ustb/data/gcnr/';   % if not found data will be downloaded from here
filename='insilico_20.uff';
tools.download(filename, url, data_path);   


mix = uff.channel_data();
mix.read([data_path filesep filename],'/mix');

%%
%% Scan
sca=uff.linear_scan('x_axis',linspace(-6e-3,6e-3,256).','z_axis', linspace(14e-3,26e-3,2.5*256).');

%%
% cyst geometry -> this should go in the uff
x0=0e-3;                
z0=20e-3; 
r=3e-3;                 

% stand off distance <- based on aperture size
M = 55;                             % aperture size
aperture = M * mix.probe.pitch;     % active aperture
F = z0 / aperture;                  % F-number
r_off = round(1.2 * mix.lambda * F, 5); % stand-off distance (rounded to 0.01 mm) 


% boundaries
ri=r-r_off;
ro=r+r_off;
rO=sqrt(ri^2+ro^2);
Ai=pi*ri^2;
Ao=pi*rO^2-pi*ro^2;
d=sqrt((sca.x-x0).^2+(sca.z-z0).^2);
%%
axi = subplot(5,4,1);
xlabel(axi,'x [mm]');ylabel(axi,'z [mm]');
viscircles(axi,[x0*1000,z0*1000],ri*1000,'EdgeColor','r','EnhanceVisibility',0);
viscircles(axi,[x0*1000,z0*1000],ro*1000,'EdgeColor','y','EnhanceVisibility',0);
viscircles(axi,[x0*1000,z0*1000],rO*1000,'EdgeColor','y','EnhanceVisibility',0);
%viscircles(axi,[x0,z0],r_speckle_outer,'EdgeColor','y','EnhanceVisibility',0);
