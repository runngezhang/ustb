% Updated 10/11/2020, Joergen Avdal (jorgen.avdal@ntnu.no)

% If using FLUST for scientific publications, please cite the original paper
% Avdal et al: Fast Flow-Line-Based Analysis of Ultrasound Spectral and
% Vector Velocity Estimators, tUFFC, 2019.

% FLUST is a simulation tool based on flowlines, useful for producing many
% realizations of the same flowfield. The motivation for making FLUST was
% to faciliate accurate assessment of statistical properties of velocity
% estimators like bias and variance.

% Using FLUST requires the user to provide/select two functions, one to create
% the flowlines, and one to calculate the PSFs from a vector of spatial
% positions. 

% After running this script, variable realTab contains realizations, and
% variables PSFs and PSFstruct contains point spread functions of the last flow line.

% FLUST accounts for interframe motion in plane wave sequences, assuming
% uniform firing rate, but not in scanned sequences.

% How to use FLUST:
% 1) Provide/select function to calculate PSFs from a vector of spatial positions. 
% 2) Run simulations with simple phantoms, check integrity of signal,
%    update quality parameters if necessary, repeat.
% 3) Run FLUST on phantom of interest.
% 4) Apply your favorite velocity estimator to realizations. 
% 5) Assess statistical properties of estimator, optimize estimator.
% 6) Publish results, report statistical properties, make results
%    reproducible.

clear all;
close all;

addpath('pathToField_II');
addpath('pathToUSTB');

s = struct();

%% DATA OUTPUT PARAMETERS
s.firing_rate = 12000; % firing rate of output signal, (Doppler PRF) = (firing rate)/(nr of firings)
s.nrReps = 10;         % nr of realizations 
s.nrSamps = 100;       % nr of slow time samples in each realization

contrastMode = 0;      % is set to 1, will simulate contrast scatterers propagating in flow field
contrastDensity = 0.1; % if using contrastMode, determines the density of scatterers, typically < 0.2

%% QUALITY PARAMETERS
s.dr = 5e-5;           % spatial discretization along flowlines: lambda/4 or smaller recommended if phase information is important
s.overSampFact = 2;    % slow time oversampling factor, should be high enough to avoid aliasing
                       % in slow time signal. Without oversampling, slow time sampling rate = firing rate

%% PERFORMANCE PARAMETER
chunksize = 5;         % chunking on scanlines, adjust according to available memory.

%% DEFINE PHANTOM AND PSF FUNCTIONS

s.phantom_function = @Phantom_small2Dtube;

s.PSF_function = @PSFfunc_L11_singlePlaneWave;
% s.PSF_function = @PSFfunc_L11_multiPlaneWave;
% s.PSF_function = @PSFfunc_LinearArray_Sectorscan;
% s.PSF_function = @PSFfunc_LinearArray_Linearscan;

s.phantom_params = []; % this structure may contain phantom parameters
s.PSF_params = [];     % this structure may contain beamforming parameters
   
%% make phantom
flowField = s.phantom_function(s.phantom_params); % flowField should have timetab and postab fields

%% FLUST main loop
for kk = 1:length(flowField)

    %% resample along flowlines with density s.dr
    prop = diff( flowField(kk).postab, 1);
    propdist = [0; cumsum( sqrt( sum( prop.^2, 2 ) ), 1 ) ];
    newdists = 0:s.dr:max(propdist);
    newtimetab = interp1( propdist, flowField(kk).timetab, newdists);
    newpostab = interp1( flowField(kk).timetab, flowField(kk).postab, newtimetab);

    %% calculate PSFs
    simStart = tic;
    PSFstruct = s.PSF_function(newpostab, s.PSF_params); % PSFs in uff/beamformed_data format
    AsimTime = toc(simStart);

    %% make realizations
    noAngs = size( PSFstruct.data, 3);
    for anglectr = 1:noAngs
        if isa( PSFstruct.scan, 'uff.sector_scan')
            szZ = length(PSFstruct.scan.depth_axis); % size( PSFs, 1);
            szX = length(PSFstruct.scan.azimuth_axis); % size( PSFs, 2);
        elseif isa( PSFstruct.scan, 'uff.linear_scan')
            szZ = length(PSFstruct.scan.z_axis); % size( PSFs, 1);
            szX = length(PSFstruct.scan.x_axis); % size( PSFs, 2);
        end

        PSFs = reshape( PSFstruct.data, [szZ, szX, noAngs, length( newtimetab)] );

        if kk == 1 && anglectr == 1
            realTab = complex( zeros( szZ, szX, s.nrSamps, noAngs, s.nrReps, 'single') );
        end

        timetab = gpuArray( newtimetab );

        ts = gpuArray( min(timetab):(1/s.firing_rate)/s.overSampFact:max(timetab) );
        Nfft = 2*length(ts)+s.nrSamps*s.overSampFact*noAngs*s.nrReps+s.overSampFact*noAngs-1;
        if anglectr == 1

            fNoiseTab = randn( [1 length( 1) length(ts)+s.nrSamps*s.overSampFact*noAngs*s.nrReps+s.overSampFact*noAngs ], 'single');

            if contrastMode
                fN_sort = sort( abs( fNoiseTab(:) ) );
                fN_thresh = fN_sort( round( length( fN_sort)*(1-contrastDensity) ) );
                fNoiseTab = single( abs(fNoiseTab) >= fN_thresh );
            end

            fNoiseTab = fNoiseTab/sqrt( length( ts) );
            fNoiseTab = fft( fNoiseTab, Nfft, 3 );

            fNoiseTab_GPU = gpuArray( fNoiseTab);
        end

        for coffset = 1:chunksize:szX
            mm = 1;
            cinds = coffset:min( coffset+chunksize-1, szX );

            myData_GPU = gpuArray( PSFs(:,cinds,anglectr,:) );

            permuteMyData = permute( myData_GPU, [4 1 2 3]);
            permuteMyData = permuteMyData(:, :);

            myData_int = interp1( timetab, permuteMyData, ts, 'linear');
            myData_int = permute( myData_int, [2 1]);
            myData_int = reshape( myData_int, [szZ, length( cinds) size( myData_int,2) ]);

            fft_myData = fft( myData_int, Nfft, 3);

            totdist = sum( sqrt( sum( diff( flowField(kk).postab,1).^2, 2 ) ), 1);
            totsamp = length( ts);
            weight = totdist/sqrt(totsamp);

            fullRealization = ifft( fft_myData .* fNoiseTab_GPU(1,1,:), [], 3);
            realTab(:,cinds,:, anglectr, : ) = realTab(:,cinds,:, anglectr, : )+...
                gather( weight*reshape( fullRealization(:,:,length(ts)+(anglectr-1)*s.overSampFact+(0:s.overSampFact*noAngs:s.nrReps*s.nrSamps*s.overSampFact*noAngs-1),:), ...
                [szZ, length( cinds), s.nrSamps, 1, s.nrReps]) );

            clc
            disp(['Flow line ' num2str(kk) '/' num2str(length( flowField) )] );
            disp(['Firing nr ' num2str(anglectr) '/' num2str(noAngs)] );
            disp(['Image line ' num2str( coffset ) '/' num2str( szX) ] );
        end
    end
end
clc
disp('Finished!')

%% VISUALIZE FIRST REALIZATION using the built inn beamformed data object
firstRealization = realTab(:,:,:,1,1);

b_data = uff.beamformed_data();
b_data.scan = PSFstruct.scan;
b_data.data = reshape(firstRealization,size(firstRealization,1)*size(firstRealization,2),1,1,size(firstRealization,3));
b_data.plot([],['Flow from FLUST'],[20])