% If using FLUST for scientific publications, please cite the original paper
% Avdal et al: Fast Flow-Line-Based Analysis of Ultrasound Spectral and
% Vector Velocity Estimators, tUFFC, 2019

% Using FLUST requires the user to provide/select two functions, one to create the
% flowlines, and one to calculate the PSFs from a vector of spatial
% positions

addpath('C:\Users\jorgenav\Documents\MATLAB\Software\Field_II');
addpath('C:\Users\jorgenav\GitProjects\ustb_flust');

%% FLUST parameters

s = struct();

s.PRF = 10000;         % PRF of output signal
s.nrReps = 10;         % nr of realizations 
s.nrSamps = 100;       % nr of slow time samples in each realization
s.overSampFact = 1;    % slow time overSampFact, high enough to avoid aliasing
s.dr = 7e-5;           % spatial discretization along flowlines: lambda/4 or smaller recommended if phase information is important

s.phantom_function = @Phantom_small2Dtube;
s.PSF_function = @PSFfunc_L11_singlePlaneWave;

s.phantom_params = []; % this structure may contain phantom parameters
s.PSF_params = [];     % this structure may contain beamforming parameters
   
%% make phantom
flowField = s.phantom_function(s.phantom_params); % flowField should have timetab and postab fields

%%
for kk = 1:length(flowField),

    %% resample along flowlines with density s.dr
    prop = diff( flowField(kk).postab, 1);
    propdist = [0; cumsum( sqrt( sum( prop.^2, 2 ) ), 1 ) ];
    newdists = 0:s.dr:max(propdist);
    newtimetab = interp1( propdist, flowField(kk).timetab, newdists);
    newpostab = interp1( flowField(kk).timetab, flowField(kk).postab, newtimetab);



    %% calculate PSFs
    simStart = tic;
    PSFs = s.PSF_function(newpostab, s.PSF_params); % PSFs in uff/beamformed_data format
    AsimTime = toc(simStart);

    %% make realizations
    anglectr = 1;

    szZ = length(PSFs.scan.z_axis); % size( PSFs, 1);
    szX = length(PSFs.scan.x_axis); % size( PSFs, 2);

    PSFs_rsh = reshape( PSFs.data, [szZ, szX, length( newtimetab)] );

    if kk == 1,
        realTab = complex( zeros( szZ, szX, s.nrSamps, size( PSFs.data, 3), s.nrReps, 'single') );
    end

    dz = PSFs.scan.z_step; %  sparams.z_axis(2)-sparams.z_axis(1);
    kz = linspace(-0.5, 0.5, szZ+1)/dz; kz(end) = [];


    incinds = [];



    chunksize = 5; %chunking on channels

    % timetab = gpuArray( flowField(pp).timetab );
    timetab = gpuArray( newtimetab );

    ts = gpuArray( min(timetab):(1/s.PRF)/s.overSampFact:max(timetab) );
    Nfft = 2*length(ts)+s.nrSamps*s.overSampFact*s.nrReps-1; %(length(ts)*(sparams.overSampFact+2) )-1;


    if anglectr == 1,
        maxoffset = 0; % do not want offset, only phase change
        offsetTab = -maxoffset:maxoffset; %linspace( -params.seed_radius, params.seed_radius, multiLineFact);


%         nrRepsNeeded = ceil( sparams.nrReps/(  length( ts)/sparams.nrSamps/sparams.overSampFact ) );
        fNoiseTab = randn( [1 length( offsetTab) length(ts)+s.nrSamps*s.overSampFact*s.nrReps ], 'single');
        fNoiseTab = fNoiseTab/sqrt( length( ts) );
        fNoiseTab = fft( fNoiseTab, Nfft, 3 );

        fNoiseTab_GPU = gpuArray( fNoiseTab);
    end


    for coffset = 1:chunksize:szX,    
        mm = 1;
        cinds = coffset:min( coffset+chunksize-1, szX );

        myData_GPU = gpuArray( PSFs_rsh(:,cinds,:) );

    %         permuteMyData = permute( myData_GPU(:, cinds, :), [3 1 2]);
        permuteMyData = permute( myData_GPU, [3 1 2]);
        permuteMyData = permuteMyData(:, :);

    %         myData_int = interp1( timetab, permuteMyData, ts, 'spline');
        myData_int = interp1( timetab, permuteMyData, ts, 'linear');
        myData_int = permute( myData_int, [2 1]);
        myData_int = reshape( myData_int, [szZ, length( cinds) size( myData_int,2) ]);


        fft_myData = fft( myData_int, Nfft, 3);


        totdist = sum( sqrt( sum( diff( flowField(kk).postab,1).^2, 2 ) ), 1);
        totsamp = length( ts);
        weight = totdist/sqrt(totsamp);
        
        fullRealization = ifft( fft_myData .* fNoiseTab_GPU(1,1,:), [], 3);
        realTab(:,cinds,:, anglectr, : ) = realTab(:,cinds,:, anglectr, : )+...
            gather( weight*reshape( fullRealization(:,:,length(ts)+(0:s.overSampFact:s.nrReps*s.nrSamps*s.overSampFact-1),:), ...
            [szZ, length( cinds), s.nrSamps, 1, s.nrReps]) );

            
        clc
        disp(['Flow line ' num2str(kk) '/' num2str(length( flowField) )] )
        disp(['Image line ' num2str( coffset ) '/' num2str( szX) ] );
    end


end