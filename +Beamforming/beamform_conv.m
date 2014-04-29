function [bf,BF] = beamform_conv(RF,varargin)

% Output:
%   bf      Beamformed and compounded signal frame
%   BF      Beamformed frames before compounding

% Rx steering is not currently implemented!!
% Spline interpolation not implemented yet.

[Nrf,Nch,Nlines,Nframes] = size(RF);

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
p.addParamValue('TxFocalDepth',10000); % planewave is default
p.addParamValue('TxAperture',[]); % in number of elements
p.addParamValue('fdemod',0);
p.addParamValue('fs',[]); % Sampling freq [Hz]
p.addParamValue('c',1540); % Sound velocity [m/s]
p.parse(varargin{:});
par = p.Results;

x = par.x;
z = par.z;
Nx = size(x,2);
Nz = length(z);

NEcorr = 10000; % make this into an input parameter
Ecorr = exp(2i*pi*single(0:NEcorr)'/NEcorr); % we need both endpoints

%% Time of arrival
% roundodd = @(v)2*round((v-1)/2)+1;
% floorodd = @(v)2*floor((v-1)/2)+1;
roundeven = @(v)2*round(v/2);
flooreven = @(v)2*floor(v/2);
% Naprt = min(floorodd(Nch), roundodd(z(:)/(par.Pitch*par.Fnumber)));
Naprt = min(flooreven(Nch), roundeven(z(:)/(par.Pitch*par.Fnumber)));

x_ele = par.x/par.Pitch + (Nch-1)/2 + 1;

%%
overshoot = max( ceil((Naprt(end)-1)/2) - (x_ele(1)-1), ...
    (x_ele(end)-1) + ceil((Naprt(end)-1)/2) - (Nch-1));
overshoot = max(0, overshoot); % make sure it's not negative

Nrf0 = 2*Nrf;
% Nchpad = [1-bfchind{end}(1),bfchind{end}(end)-Nch];
RF0 = zeros(Nrf0,Nch+2*overshoot,Nlines,'like',RF); % make zero-padded frame
RF0_ = RF0;
bf = zeros(Nz,Nx,'like',RF);

%% Apodization and Rx steering
[Naprt,~,uind] = unique(Naprt);
for kk = 1:numel(Naprt),
    nn = -(Naprt(kk)-1)/2:(Naprt(kk)-1)/2;
    bfch = bsxfun(@plus,nn(:),overshoot+x_ele);
    % channels used in beamforming each line
    bfch = fix_rounding_problems(bfch,1e-6);
    bfchind{kk} = bsxfun(@plus, round(bfch), (0:Nlines-1)*size(RF0,2));
    apod{kk} = generateApodization(par.Apodization,Naprt(kk));
end

Lhalf = (par.TxAperture-1)/2 * par.Pitch;
tau_foc = (sqrt(par.TxFocalDepth^2 + Lhalf.^2) - par.TxFocalDepth) / par.c;
for kk = 1:Nz,
    Ahalf = (-(Naprt(uind(kk))-1)/2:(Naprt(uind(kk))-1)/2) * par.Pitch;
    toa{kk} = z(kk)/par.c + tau_foc + 1/par.c * sqrt(z(kk)^2 + Ahalf(:).^2) - par.RxDelay;
    tind{kk} = bsxfun(@plus, (bfchind{uind(kk)}-1)*Nrf0, toa{kk}*par.fs);
    phase_corr{kk} = par.fdemod*(toa{kk}-2*z(kk)/par.c);
end



%% This is the main loop
for jj = 1:Nframes,
    RF0(1:Nrf,overshoot+(1:Nch),:) = RF(:,:,:,jj);
    RF0_(1:Nrf-1,overshoot+(1:Nch),:) = diff(RF(:,:,:,jj));
    for kk = 1:Nz,
        awin = repmat(apod{uind(kk)},1,ceil(Nx/size(apod{uind(kk)},2)));
        bf(kk,:,jj) = sum( awin(:,1:Nx) .* ...
            interp_local(single(tind{kk}),RF0,RF0_,phase_corr{kk},par.Interpolation));
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
            yi = yi .* repmat(Ecorr(eind),1,size(yi,2)); % move this out
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