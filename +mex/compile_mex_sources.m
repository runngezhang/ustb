% compile mex sources for

% linux
%mex source/cpwir.cpp -I"/media/Data/ingvilek/ustb2/Source Code/CPWI_LR/Mex/tbb44_20151115oss/include" -L"/media/Data/ingvilek/ustb2/Source Code/CPWI_LR/Mex/tbb44_20151115oss/lib/intel64/gcc4.4" -ltbb

copyfile(fullfile(matlabroot,'extern'),'.mex','f')

%% CPWC 
disp('------------------------ CPWC');
mex source/cpwir.cpp

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
mex source/snell.cpp
mex source/snell_detector.cpp