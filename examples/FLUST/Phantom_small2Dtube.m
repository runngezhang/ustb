function flowField = Phantom_small2Dtube( p ) % parameter structure p not used in this example

%% small 2D tube phantom for quick demo purposes
btf = 60;
npoints = 10;
timeend = 0.01; %0.0075;
timecent = timeend/2;
tubedepth = 0.02;
step = 0.00015; %lambda/2 for 5 MHz
max_vel = 1;
for kk = 1:20,
    currtubedepth = tubedepth + step*(kk-1);
    flowField(kk).timetab = linspace(0, timeend, npoints);
    flowField(kk).postab = max_vel*(flowField(kk).timetab-timecent).*[sind(btf); 0; cosd(btf)]+[0; 0; currtubedepth];
    flowField(kk).timetab = flowField(kk).timetab.'; 
    flowField(kk).postab = flowField(kk).postab.';
end