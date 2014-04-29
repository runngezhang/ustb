function s = setup_tci()
    %%
    s.basic.c    = 1540;
    s.basic.framenumber = 1;
    s.basic.fc   =   containers.Map({'Viglen',       'Ruten',       'Vora'},...
                                    {[6]*1e6,        [7]*1e6,      [8]*1e6});       
    
    s.basic.fileinfo_review = false;
    s.basic.fileinfo_showall = false;
    s.basic.fileinfo_params = [IdDb.FOCUS_DEPTH, IdDb.IQ_SAMPLE_FREQ, IdDb.FREQ, IdDb.LF_FREQ_SURF, IdDb.LF_VOLTAGE_SURF, IdDb.LF_HALF_CYCLES_SURF, ...
                               IdDb.LF_INDEP_FOCUS_SURF, IdDb.LF_FOC_DEPTH_SURF, IdDb.LF_CENTER_ACTIVE_SURF, IdDb.LF_OUTER_ACTIVE_SURF, IdDb.HF_POSITION_SURF...
                               1334:1344];
    s.basic.fileformat = '0.1';  % old ultrasonix file format is default
    s.basic.setupname = 'setup_tci';
    
    % RF filter
    s.rffilter.interlace = true; % drop every second line?
    s.rffilter.order = 52;
    s.rffilter.cutoff = containers.Map({'Viglen',       'Ruten',        'Vora'},...
                                       {[3.5 9.5]*1e6,  [3  12]*1e6,   [3.4  14]*1e6});
    s.rffilter.review = false;
    
  %%%%%%%%%%%%%%%%
  %IQ demodulation                                   
    % B-Mode demodulation frequency
    s.iqdemod.fdemodB = containers.Map({'Viglen',       'Ruten',      'Vora'},...
                                       {8*1e6,        8*1e6,           7e6});
    % B-Mode demodulation cutoff frequency for narrowband demod
    s.iqdemod.cutoffB = containers.Map({'Viglen',       'Ruten',      'Vora'},...
                                       {2.5*1e6,          3*1e6,      2.2e6});
    % SURF-Mode demodulation frequency
    s.iqdemod.fdemodS = containers.Map({'Viglen',       'Ruten',      'Vora'},...
                                       {8*1e6,        8*1e6,          7e6});
    % SURF-Mode demodulation cutoff frequency for narrowband demod
    s.iqdemod.cutoffS = containers.Map({'Viglen',       'Ruten',      'Vora'},...
                                       {2.5*1e6,          3*1e6,      2.2e6});
    % filter order for narrow band demod
    s.iqdemod.order = 35;
    s.iqdemod.iq_dec = 1;   % decimation, 1= no decimation, 2: every other sample is left out
    s.iqdemod.review = false;
%     s.iqdemod.review = true;
    
  %%%%%%%%%%%%%%%%%%%%%  
  % delay estimation
    s.delayEst.win_mm = .4;
    s.delayEst.reg_meth = 'const_conv';
    s.delayEst.filterBilat = true;
%     s.delayEst.filterBilat = false;
    
  %%%%%%%%%%%%%%%%%%%%%
  % delay correction
    s.delayCorr.tau_factor = 0.50;  % tau/half by default
    s.delayCorr.tau_iteration_start = 1;
    s.delayCorr.tau_thres_fact = 1; % gain factor for threshold (compensated for tau_gain, so 1 gives same as the correction delay)
%     s.delayCorr.tau_thres_fact = 1; % gain factor for threshold (compensated for tau_gain, so 1 gives same as the correction delay)
   %   s.delayCorr.tau_iteration_start = 1;
    s.delayCorr.tau_gain = 0.9; % gain to adjust global level of delay which is used for correction
%     s.delayCorr.tau_gain = 1; % gain to adjust global level of delay which is used for correction
    s.delayCorr.method = 'wall';
%     s.delayCorr.method = 'no_avg';
%     s.delayCorr.method = 'no_avg_to0';
%     s.delayCorr.method = 'lin_extrap';
%       s.delayCorr.method = 'lin_extrap_first200';
      s.delayCorr.method = 'lin_extrap_first200to0';
%       s.delayCorr.method = 'lin_extrap_first200to0_llapprox';
%       s.delayCorr.method = 'lin_extrap_first200to0_llapprox_lat';
%       s.delayCorr.method = 'lin_extrap_first200to0-5lines';
%       s.delayCorr.method = 'minall';
%     s.delayCorr.method = 'apriori-phantom-file';
    
    s.delayCorr.correctProbeDefects = 'no';
%     s.delayCorr.correctProbeDefects = 'Viglen_0';
%     s.delayCorr.correctProbeDefects = 'Vora_0';
    
    
    s.delayCorr.file = '';
    s.delayCorr.loadMean = []; % just for apriori-phantom-file 
    s.delayCorr.interpMethod = 'linear';
    
    s.delayCorr.review = true; % show plots on estimated delay and delay correction
    s.delayCorr.review = false;
  %%%%%%%%%%%%%%%%%%%%%
  % tci processing methods
<<<<<<< .mine
    s.tciProc.method = 'phaseweight'; % one of subtract, tau/2, phaseweight or phasedifference                             
    s.tciProc.method = 'plus';
    s.tciProc.method = 'tau_half';
=======
%     s.tciProc.method = 'phase_weight'; % one of subtract, tau/2, phaseweight or phasedifference
                             %see rev_suppression.m for details
%      s.tciProc.method = 'plus-iter';
     s.tciProc.method = 'plus';
%      s.tciProc.method = 'plus-noavg';
%      s.tciProc.method = 'plus-test';
>>>>>>> .r555
%      s.tciProc.method = 'plus-integral';
%     s.tciProc.method = 'plus-thres';
%     s.tciProc.method = 'tau_half';
%     s.tciProc.method = 'tau_x';
%     s.tciProc.method = 'tau_half_integral';
    s.tciProc.win_m = 0.5*1e-3; % window size for delay median filter if applicable
    s.tciProc.interpMethod = 'linear'; % interpolation method if delay processing is done once more    
    
  %%%%%%%%%%%%%%%%%%%%%
  % lateral upsampling 
    s.latUpsample.order = 2;
    s.latUpsample.dim = 2;
    
  %%%%%%%%%%%%%%%%%%%%%
  % Gain compensation
    s.gainComp.enable = true;   % do gain compensation?
%     s.gainComp.enable = false;   % do gain compensation?
    switch s.tciProc.method
        case 'tau_half'
           s.gainComp.method = 'parametric-filt-half';            
        case 'plus'
           s.gainComp.method = 'parametric-filt-plus';
        case 'phase_weight'
           s.gainComp.method = 'parametric-filt-plus';
        case 'plus-thres'
           s.gainComp.method = 'parametric-filt-plus';
        case 'plus-test'           
           s.gainComp.enable = false;   % do gain compensation?
        case 'plus-integral'
           s.gainComp.method = 'parametric-filt-plus';
        otherwise
           s.gainComp.method = 'latavg';
    end;
%     s.gainComp.method = 'latavg';
%     s.gainComp.method = 'parametric';
%     s.gainComp.method = 'parametric-filt-plus';
%     s.gainComp.method = 'parametric-filt-half';
    s.gainComp.reviewGain = -20;
    s.gainComp.reviewDyn = 60;
    s.gainComp.review = true;
    s.gainComp.review = false;

    
  %%%%%%%%%%%%%%%%%%%%%
  % scan conversion
    s.scanConv.Nx = 512;
    s.scanConv.Nz = 512;
    s.scanConv.keepAR = true; % if true, Nz is ignored and chosen according to the aspect ratio of the image
    s.scanConv.review = false;
    s.scanConv.review = true;
    
  %%%%%%%%%%%%%%%%%%%%%
  % image compression
    
    % default gain

    s.image.gainB = -38;
    s.image.dynB = 60;
    s.image.gainS = -38;
    s.image.dynS = 60;

    s.image.max = 255;
    s.image.value_max = 2^16-1;
    s.image.review = true;
    s.image.review = false;

    s.image.doSpeckleReduction = true;
    %s.image.doSpeckleReduction = false;
    
    s.image.gcurveB = [];
    s.image.gcurveS = [];
    s.image.gainCurveFile = '';    
    s.image.takeGainCurveSnapshot = true;
    s.image.gainCurveWriteBack = true; % if the gainfile is present, save the gain curve if changed?
    s.image.gainCurveWriteBack = false; % if the gainfile is present, save the gain curve if changed?
    
    %switch parent.tciProc.method
    %    case 'tau_half'
    %    case 'phase_weight'
    %    case 'plus'
    %    case 'plus-thres'    
    %end;
    
  %%%%%%%%%%%%%%%%%%%%%
  % persistence
    s.persist = struct();
    s.persist.enable = true;
    s.persist.length = 3;
    s.persist.window = 'hamming';
    
  %%%%%%%%%%%%%%%%%%%%%
  % movie
    s.movie = struct();
    s.movie.enable = true;
    s.movie.enable = false;
    s.movie.batchMode = false; % don't start gaining
    s.movie.batchpreview = false; % just process 1 image for preview
    s.movie.batchnogaindyn = false;
    s.movie.batchnogaindyn = false;
    s.movie.frames = [1:20];
    s.movie.fps = [];
    s.movie.postfix = '';
    