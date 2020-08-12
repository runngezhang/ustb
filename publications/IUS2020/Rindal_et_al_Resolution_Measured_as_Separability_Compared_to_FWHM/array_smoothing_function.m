filename='L7_FI_IUS2018.uff';
FI_channel_data = uff.read_object([data_path filesep filename],'/channel_data');

%%

FI_channel_data.probe.pitch

probe = FI_channel_data.probe;





%%
%clear all; close all
%Defining lambda
%p = parameters(0,1,27,128);
%lambda = probe.lambda;
d = probe.pitch;
D = probe.N_elements * d;
elemt_size = probe.width(1);
 
 
%Creating Kx wector
kx = linspace(-3*pi/d,3*pi/d,1000);
kx_center =linspace(-(0.5)/d,(0.5)/d,2000);
for i = 1:probe.N_elements
    x(i) = i*d - (probe.N_elements/2)*d;
end
w_m = ones(1,probe.N_elements);
% Calculate the Array pattern
% Input :
%       kx  : wavenumber
%       x   : element position
%       w_m : element weighting
 
W = zeros(length(kx),1);
for t = 1:length(kx)
    W(t) = sum(w_m.*exp(-1i*kx(t).*x));
end
 
W_center = zeros(length(kx),1);
for t = 1:length(kx_center)
    W_center(t) = sum(w_m.*exp(-1i*kx_center(t).*x));
end
 
W_el =  sin(kx.*(elemt_size/2))./(kx/2);
W_total = W_el.*W';
 
W_el_center = sin(kx_center.*(elemt_size/2))./(kx_center/2);
W_total_center = W_el_center.*W_center';
 
visible= 2*pi/lambda;
 
h1 = figure(1)
clf
%subplot(211)
% plot(x*1000,zeros(1,length(x)),'*')
% xlabel('mm')
% title('Element positions');
 
subplot(211)
plot(kx*d,db(W)-max(db(W)));
hold all
plot(kx*d,db(W_total)-max(db(W_total)))
plot(kx*d,db(W_el)-max(db(W_el)));
plot([-visible*d, -visible*d],[-50 0],'-r','LineWidth',2);
plot([visible*d, visible*d],[-50 0],'-r','LineWidth',2);
xlabel('Wavenumber k_xd');
axis([-3*pi 3*pi -50 0]);
title('Aperture function');
ylabel('db(W)');
legend('Warray','Welement','Wtot');
 
 
subplot(212)
dbW = db(W_center);
dbNorm = dbW-max(dbW(:));
%plot(kx_center*d,dbNorm);
hold all
plot(kx_center*d,db(W_total_center)-max(db(W_total_center)),'LineWidth',2)
plot((kx_center*d),-3*ones(1,length(kx_center)),'LineWidth',2)
plot(kx_center*lambda,-6*ones(1,length(kx_center)),'LineWidth',2)
%plot(kx_center*d,db(W_el_center)-max(db(W_el_center)));
xlabel('Wavenumber k_xd');
title('Aperture function zoomed')
ylabel('db(W)');
axis([-0.2 0.2 -20 0]);
legend('Wtotal','-3 db','-6 db');
set(findall(h1,'type','text'),'fontSize',18)
set(gca(h1),'fontSize',14)

angular_x_axis = D.*asin(lambda.*kx_center/(2*pi));

figure(101);clf;
subplot(211);hold on;
plot(angular_x_axis,db(W_total_center)-max(db(W_total_center)));
plot(angular_x_axis,db(W_center)-max(db(W_center)),'LineWidth',2);
plot((angular_x_axis),-3*ones(1,length(kx_center)),'LineWidth',2)
plot(angular_x_axis,-6*ones(1,length(kx_center)),'LineWidth',2)
res = ((1.22*lambda*D)/(D))
hold all;
plot([-res/2 -res/2],[-60 0])
plot([res/2 res/2],[-60 0])
plot(scan.x_axis,img_field_II(69,:)); 

subplot(212)
plot(db(W_total_center)-max(db(W_total_center)));

%%
db_W_two_way = db(W_total_center) - max(db(W_total_center(:)));
[V,I] = sort(abs((db_W_two_way + 3)),'ascend');
two_way_3dbres = abs(kx_center(I(1)))
[V,I] = sort(abs((db_W_two_way + 3/2)),'ascend');
two_way_15dbres = abs(kx_center(I(1)))
[V,I] = sort(abs((db_W_two_way + 3)),'ascend');

%% v

%value = 0.02226
value = 0.02977;
angular_test = rad2deg(lambda.*0.0215/(2*pi))
angular = rad2deg(asin(lambda.*value/(2*pi)))

angular.*(probe.N*probe.pitch)

D*(1.22*lambda/D)*10^3

1.22*lambda*10^3

tan(angular/2)*D


%%

angle = lambda/D;

res = lambda*D/D*10^3