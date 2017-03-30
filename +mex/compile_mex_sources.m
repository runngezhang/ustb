% compile mex sources for

% Works on MS with VS2010+ and Linux with Intel tbb

copyfile(fullfile(matlabroot,'extern'),'.mex','f')

%% CPWC 
disp('------------------------ CPWC');
mex.build_mex('cpwir');

%% CPWC Low resolution images
disp('------------------------ Low res CPWC');
mex.build_mex('cpwlr');

%% CPWC Low-Low resolution images
disp('------------------------ Low-low res CPWC');
mex.build_mex('cpwllr');

%% compile thor-andreas code
disp('------------------------ Thor-Andreas CPWC');
mex.build_bf(1,1,1,1);

%% compile snell
disp('------------------------ Snell beamformer STAI');
mex.build_mex('snell');

%% compile snell-detector
disp('------------------------ Snell detector');
mex.build_mex('snell_detector');
