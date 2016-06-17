% compile mex sources for

% linux
%mex source/cpwir.cpp -I"/media/Data/ingvilek/ustb2/Source Code/CPWI_LR/Mex/tbb44_20151115oss/include" -L"/media/Data/ingvilek/ustb2/Source Code/CPWI_LR/Mex/tbb44_20151115oss/lib/intel64/gcc4.4" -ltbb

copyfile(fullfile(matlabroot,'extern'),'.mex','f')

%% CPWC 
mex source/cpwir.cpp

%% CPWC Low resolution images
mex.build_mex('cpwlr');

%% CPWC Low-Low resolution images
mex.build_mex('cpwllr');

%% compile thor-andreas code
mex.build_bf(1,1,1,1)