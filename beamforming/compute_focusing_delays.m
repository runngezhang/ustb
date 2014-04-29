function tau = compute_focusing_delays(z,x0,c0,th)

z = z(ones(1,length(x0)),:)';
x0 = x0(ones(1,length(z)),:);

xF = z*sin(deg2rad(th));
yF = z*cos(deg2rad(th));

tau = (sqrt((xF-x0).^2 + yF.^2) - z)/c0;    