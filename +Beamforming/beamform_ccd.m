function [bf,BF] = beamform_ccd(RF,varargin)

% Output:
%   bf      Beamformed and compounded signal frame
%   BF      Beamformed frames before compounding

% Rx steering is not currently implemented!!
% Spline interpolation not implemented yet.

[Nrf,Nch,Nframes] = size(RF);

%%
p = inputParser();
p.addOptional('x',-(Nch-1)/2:(Nch-1)/2); % spatial positions azimuth [m]
p.addOptional('z',[]); % spatial sampling position depth [m]
p.addOptional('TxSteeringAngle',0); % [rad]
p.addOptional('RxSteeringAngle',0); % use this to modify x and y
% p.addOptional('RxDepthOffset',0); % use this to modify x and y
p.addOptional('RxDelay',0); % Delay before beginning to sample [s]
p.addParamValue('Fnumber',0); % receive f-number
p.addParamValue('Apodization',@hamming); % only hamming accepted right now
p.addParamValue('Interpolation','linear');
p.addParamValue('Pitch',[]);
p.addParamValue('fdemod',0);
p.addParamValue('fs',[]); % Sampling freq [Hz]
p.addParamValue('c',1540); % Sound velocity [m/s]
p.parse(varargin{:});
par = p.Results;

NEcorr = 10000; % make this into an input parameter
Ecorr = exp(2i*pi*single(0:NEcorr)'/NEcorr); % we need both endpoints

%% Not taking account of Rx steering yet
x0 = (-(Nch-1)/2:(Nch-1)/2) * par.Pitch;
xmin = x0(1);

x = par.x;
Nx = size(x,2);
xind = 1:numel(x0); % 1 + round((x-x0(1))/par.Pitch); % round this?
z = par.z(:);
Nz = numel(z);

% projection of each grid point onto transducer (xxd) along steering angle
phi = par.RxSteeringAngle;
if phi ~= 0,
    xxd = bsxfun(@minus,x,z*tan(phi));
else
    xxd = repmat(x,length(z),1);
end
xxd_ele = xxd/par.Pitch + (Nch+1)/2;

%% Apodization and Rx steering
Naprt = min(Nch, round(z(:)/(par.Pitch*par.Fnumber)));
overshoot = ceil( Naprt(end)/2 + max(-min(xxd_ele(:)), max(xxd_ele(:))-Nch));
overshoot = max(0, overshoot); % avoid negative values
x0padded = [(-overshoot:-1)*par.Pitch+x0(1),x0(:)',(1:overshoot)*par.Pitch+x0(end)];
for kk = 1:numel(Naprt),
    nn = -(Naprt(kk)-1)/2:(Naprt(kk)-1)/2;
    bfch = bsxfun(@plus,nn(:),overshoot+xxd_ele(kk,:));
    % channels used in beamforming each line
    bfch = fix_rounding_problems(bfch,1e-6);
    bfchind{kk} = round(bfch);
    bfch_shift = bfchind{kk}(1,:) - bfch(1,:);
    % apod{kk} = generateApodization(par.Apodization,Naprt(kk));
    % xxd_shift_unique = almostunique(bfch_shift,1e-4);
    xxd_shift_unique = unique(bfch_shift);
    % assert(numel(xxd_shift_unique)==1,'Doesn''t yet work for these x-values.');
    for ii = 1:numel(xxd_shift_unique),
        apod{kk}(:,ii) = generateApodization(par.Apodization,Naprt(kk),xxd_shift_unique(ii));
    end
end

%% 
% this will really only work if x is the same as the positions of all
% the channels. Should fix this in case we want e.g. narrower sector
Nrf0 = 2*Nrf;
% Nchpad = [1-bfchind{end}(1),bfchind{end}(end)-Nch];
RF0 = zeros(Nrf0,Nch+2*overshoot,'like',RF); % make zero-padded frame
RF0_ = RF0;
bf = zeros(Nz,Nx,'like',RF);

return_all_frames = (nargout > 1);
% Also return beamformed data from individual plane wave
if return_all_frames,
    BF = zeros(Nz,Nx,Nframes,'like',RF);
end

d2 = cell(Nz,1);
for mm = 1:Nz,
    % some of this is redundant without Rx steering, but takes very little
    % time.
    x_ = x(min(mm,end),:);
    x2 = bsxfun(@minus,x0padded(:),x_).^2;
    x2ind = bsxfun(@plus,bfchind{mm},(0:Nx-1)*(Nch+2*overshoot));
    d2{mm} = sqrt(z(mm)^2 + x2(x2ind));
end

%% This is the main loop
for jj = 1:Nframes,
    RF0(1:Nrf,overshoot+(1:Nch)) = RF(:,:,jj);
    RF0_(1:Nrf-1,overshoot+(1:Nch)) = diff(RF(:,:,jj));
    theta = par.TxSteeringAngle(mod(jj-1,end)+1);
    for kk = 1:Nz,
        x_ = x(min(kk,end),:);
        d1 = (x_*sign(theta)-xmin)*sin(abs(theta)) + z(kk)*cos(theta);
        toa = bsxfun(@plus,d1/par.c,d2{kk}/par.c) - par.RxDelay;
        tind = (bfchind{kk}-1)*Nrf0 + toa*par.fs;
        phase_corr = par.fdemod*(toa-2*z(kk)/par.c);
        % bf0 = apod{kk}' * ...
        %     interp_local(single(tind),RF0,RF0_,phase_corr,par.Interpolation);
        awin = repmat(apod{kk},1,ceil(Nx/size(apod{kk},2)));
        bf0 = sum( awin(:,1:Nx) .* ...
            interp_local(single(tind),RF0,RF0_,phase_corr,par.Interpolation));
        % assert(all(almostequal(bf0(:),bf0_(:),1e-3)),'Uh-oh!');
        bf(kk,:) = bf(kk,:) + bf0;
        if return_all_frames,
            BF(kk,:,jj) = bf0;
        end
    end
end

%% Nested functions

    function yi = interp_local(t,y0,y0diff,phi,method)
        switch lower(method),
            case 'nearest'
                yi = y0(round(t));
            case 'linear',
                t0 = floor(t);
                dt = t - t0;
                t0int = uint32(t0);
                % yi = y0(t0int).*(1-dt) + y0(t0int+1).*dt;
                % yi = y0(t0int) + (y0(t0int+1)-y0(t0int)).*dt;
                yi = y0(t0int) + y0diff(t0int).*dt;
                % yi = yi .* exp(2i*f_demod*(t_ix - dec*(iz-1)))
            case 'spline',
                error('Spline interpolation is currently not supported');
            otherwise
                error(['Unknown interpolation type: %s.\n', ...
                    'Use ''nearest'' or ''linear''.'],method);
        end
        if par.fdemod,
            %E = exp(2i*pi*phi);
            % yi = yi .* E;
            % Using an approximation here with table lookup
            eind = uint32((phi-floor(phi))*NEcorr)+1;
            yi = yi .* Ecorr(eind);
        end
    end

%% 
    function w = generateApodization(wintype,winlength,wshift)
        if nargin > 2,
            for qq = 1:numel(wshift),
                w = myhamming(winlength,wshift(qq));
            end
            return
        end
        try
            w = window(wintype,winlength);
        catch err,
            if strcmpi(err.identifier,'MATLAB:UndefinedFunction'),
                cstr = regexp(help('window'),'\n','split');
                cmatch = regexp(cstr,'(?<=^\s*)@\w+','match','once');
                error(['Unknown apodization window: "%s".\n', ...
                    'Please select one of the following window types:\n%s'], ...
                    ['@',char(wintype)],[char(cmatch),repmat(' ',numel(cmatch),1)]');
            else
                err.rethrow();
            end
        end
    end

%% Approximate array with rational values to avoid floating point round-off errors
    function value_nice = fix_rounding_problems(value_naughty,varargin)
        % Approximate floating point value with 'nice' rational numbers, to
        % make rounding more pleasant. Hoping to avoid numbers differing by
        % ~eps.
        
        % Much too slow!!
        % [N,D] = rat(value_naughty,varargin{:});
        % value_nice = N./D;
        
        % almostunique fails on the return of JJ right now. Consider a fix.
        % [vnu,~,JJ] = almostunique(value_naughty(:),varargin{:});
        vn_floor = floor(value_naughty(:));
        vn_rem = value_naughty(:) - vn_floor;
        [vn_rem_u,~,JJ] = unique(vn_rem); % a bit slower in some cases
        [N,D] = rat(vn_rem_u,varargin{:});
        value_nice = reshape(vn_floor + N(JJ)./D(JJ), size(value_naughty));
    end


end