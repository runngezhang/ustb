    s.basic = struct();
    s.rffilter = struct();
    s.iqdemod = struct();
    s.delayEst = struct();
    s.delayCorr = struct();
    s.tciProc = struct();
    s.latUpsample = struct();
    s.gainCompIQ = struct();
    s.scanConv = struct();
    s.image = struct();
    s.persist = struct();
    s.movie = struct();
    
%  NB: The values in the comments here may have changed from the global
%  file and are just given for help how to define the different parameters


%   %%%%%%%%%%%%%%%%%%%%%
%   % image compression
% 
%     % default gain
<<<<<<< .mine
%    switch s.tciProc.method
%        case 'plus'
%        case 'phaseweight'
%        case 'tau_half'
%    end;
=======
    switch parent.tciProc.method
        case 'tau_half'
        case 'phase_weight'
        case 'plus'
        case 'plus-thres'
    end;
>>>>>>> .r555
%     s.image.gainB = -20;
%     s.image.dynB = 60;
%     s.image.gainS = -20;
%     s.image.dynS = 60;
%     s.image.max = 255;
%     s.image.value_max = 2^16-1;
%     s.image.review = true;


%     s.basic.c    = 1540;
%     s.basic.framenumber = 10;
%     s.basic.fc   =   containers.Map({'Viglen',       'Ruten'},...
%                                     {[6]*1e6,        [7]*1e6});       
% 
%     s.basic.fileinfo_review = true;
%     s.basic.fileinfo_showall = false;
%     s.basic.fileinfo_params = [IdDb.FOCUS_DEPTH, IdDb.IQ_SAMPLE_FREQ, IdDb.FREQ, IdDb.LF_FREQ_SURF, IdDb.LF_VOLTAGE_SURF, IdDb.LF_HALF_CYCLES_SURF, ...
%                                IdDb.LF_INDEP_FOCUS_SURF, IdDb.LF_FOC_DEPTH_SURF, IdDb.LF_CENTER_ACTIVE_SURF, IdDb.LF_OUTER_ACTIVE_SURF, IdDb.HF_POSITION_SURF...
%                                1334:1344];
%     % RF filter
%     s.rffilter.order = 52;
%     s.rffilter.cutoff = containers.Map({'Viglen',       'Ruten'},...
%                                        {[3.5 9.5]*1e6,  [3  12]*1e6});
%     s.rffilter.review = true;
% 
%   %%%%%%%%%%%%%%%%
%   %IQ demodulation                                   
%     % B-Mode demodulation frequency
%     s.iqdemod.fdemodB = containers.Map({'Viglen',       'Ruten'},...
%                                        {8*1e6,        8*1e6});
%     % B-Mode demodulation cutoff frequency for narrowband demod
%     s.iqdemod.cutoffB = containers.Map({'Viglen',       'Ruten'},...
%                                        {2.5*1e6,          3*1e6});
%     % SURF-Mode demodulation frequency
%     s.iqdemod.fdemodS = containers.Map({'Viglen',       'Ruten'},...
%                                        {8*1e6,        8*1e6});
%     % SURF-Mode demodulation cutoff frequency for narrowband demod
%     s.iqdemod.cutoffS = containers.Map({'Viglen',       'Ruten'},...
%                                        {2.5*1e6,          3*1e6});
%     % filter order for narrow band demod
%     s.iqdemod.order = 50;
%     s.iqdemod.iq_dec = 1;   % decimation, 1= no decimation, 2: every other sample is left out
%     s.iqdemod.review = true;
% 
%   %%%%%%%%%%%%%%%%%%%%%  
%   % delay estimation
%     s.delayEst.win_mm = 1;
%     s.delayEst.reg_meth = 'rlin_conv';
% 
%   %%%%%%%%%%%%%%%%%%%%%
%   % delay correction
%     s.delayCorr.tau_factor = 0.5;  % tau/half by default
%     s.delayCorr.tau_gain = 1; % gain to adjust global level of delay which is used for correction
%     s.delayCorr.method = 'wallViglen';
% %     s.delayCorr.method = 'no_avg';
% 
% %     s.delayCorr.method = 'apriori-phantom-file';
%     s.delayCorr.file = '';
%     s.delayCorr.loadMean = []; % just for apriori-phantom-file 
% 
%     s.delayCorr.review = true; % show plots on estimated delay and delay correction
% 
%   %%%%%%%%%%%%%%%%%%%%%
%   % lateral upsampling 
%     s.latUpsample.order = 1;
%     s.latUpsample.dim = 2;
% 
%   %%%%%%%%%%%%%%%%%%%%%
%   % Gain compensation
%     s.gainComp.enable = true;   % do gain compensation?
%     s.gainComp.reviewGain = -20;
%     s.gainComp.reviewDyn = 60;
%     s.gainComp.review = true;
%   %%%%%%%%%%%%%%%%%%%%%
%   % scan conversion
% 
%     s.scanConv.Nx = 640;
%     s.scanConv.Nz = 480;
%     s.scanConv.review = true;
% 
% 
%   %%%%%%%%%%%%%%%%%%%%%
%   % persistence
%     s.persist = struct();
%     s.persist.enable = true;
%     s.persist.length = 7;
%     s.persist.window = 'hamming';
% 
%   %%%%%%%%%%%%%%%%%%%%%
%   % movie
%     s.movie = struct();
%     s.movie.enable = true;
% %     s.movie.enable = false;
%     s.movie.frames = [];
%     s.movie.fps = [];
% s.movie.postfix = '';