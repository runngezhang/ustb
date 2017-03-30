% compile mex sources for

% FOR WINDOWS:
% Works on MS with Visual Studo 2010+.

% FOR LINUX :
% All code exept Thor Andreas Code compiles if you have the Intel tbb
% installed

% FOR MAC:
% All code exept Thor Andreas Code compiles if you have installed the 
% Intel TBB library. This is easiest done with 
% "brew install tbb". Thus, using the Homebrew (https://brew.sh/index_no.html). 
% It is then defaultly installed to /usr/local/Cellar/tbb/2017_U5/include/tbb/

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
