function [bf X Z t_ix ind apod] = planewave_beamforming(ch_data,p,x,z,type,apod_type,varargin)
% bf = planewave_beamforming(ch_data,p,x,z,type,apod_func)
% 
% Beamforming method for plane wave emissions. The parameter structure 'p'
% should have the follwoing fields:
% c         - Speed f sound [m/s]
% pitch     - Array pitch
% kerf      - Array kerf
% t0        - Time of first sample
% tx_angle  - Transmit angle in radians
% rx_angle  - Receive angle in radians
% FN        - F number for receive aperture.
% fs_in     - Sampling frequency of input signal
% fs_ou     - Sampling frequency of output
% f_demod   - IQ demodulation frequency. Needed if data is complex
% 
% Input:
% rf_ch     - Channel data, RF or IQ
% p         - Parameter structure
% x         - x coordinates for beamforming area. If receive steering is
% applied this is the x coordinate of the starting point if each receive
% line. The offset is computed
% z         - z coordinates for beamforming area. If receive steering is
% applied this is the radial coordinates.
% type      - Interpolation type, 'nearest', 'linear' or 'cubic'
% apod_func - [Optional] Apodization. Default: 'rect'. Other values
% 'hamming', 'tukey', 'cheb'. For tukey and cheb, also specify window
% parameter
%
% Output
% bf        - Beamformed data.
% X         - x-grid coordinates for the beamformed data
% Z         - z-grid coordinates for the beamformed data
% 

% Set default apodization to rectangular
if nargin < 6 || isempty(apod_type)
    apod_type = 'rect';
end

% Depending on if data is real or not, do some stuff
if isreal(ch_data) 
    if isfield(p,'f_demod') && p.f_demod > 0
        warning('If input data is real, f_demod should be zero. Setting f_demod to zero.')
    end
    p.f_demod = 0;      
end

if isfield(p,'fs') || isfield(p,'fs_iq')
    error('fs and fs_iq is deprecated. Use fs_in and fs_out. See documentation');
end

% Input data dimensions
p.nframes   = size(ch_data,3);
p.nch       = size(ch_data,2);
p.nsamples  = size(ch_data,1);

% x coordinate of transducer elements
x_ch    = (0:(p.nch-1))'*p.pitch;

if abs(p.rx_angle) > 0
    dx      = z'*sin(p.rx_angle);
    z       = z*cos(p.rx_angle);
    
    [X0 Z]  = meshgrid(x,z);
    X       = X0 + dx(:,ones(1,length(x)));

    X       = X(:,:,ones(1,p.nch));
    Z       = Z(:,:,ones(1,p.nch));
    x_ch    = reshape(x_ch,1,1,[]);
    X_ch    = x_ch(ones(1,size(X,1)),ones(1,size(X,2)),:);
    
    % Compute aperture. Aperture is choosen based on the starting
    % x-coordinate of the beamforming line.    
    % apod = abs((X_ch-X0(ones(1,length(z)),:,ones(1,p.nch)))./Z)/(1./p.FN/2);    
else
    [X Z X_ch] = meshgrid(x,z,x_ch);            
    
    % apod = (abs(Z./(X_ch-X)/2))/(p.FN);    
end

% Compute apodization
%apod = gen_apodization(apod,apod_type,varargin{:});

apod = generate_apodization(p,z,x-p.pitch*(p.nch-1)/2,p.c/p.fc,apod_type);

Nz = 201;
Nx = 3;
h = ones(Nz,Nx)/Nx/Nz;
apod = single(imfilter(apod,h));

% Dimensions and matrix allocation
nx  = size(X,2);
nz  = size(X,1);
bf  = zeros(nz,nx,p.nframes);                  

% Depth and channel indices
iz      = (1:nz)';
iz      = iz(:,ones(1,nx),ones(1,p.nch));
ch_ix   = reshape((1:p.nch)',1,1,[]);
ch_ix   = ch_ix(ones(1,nz),ones(1,nx),:);

% Angular demodulation frequency, normalized to sampling frequency
w_demod = 2*pi*p.f_demod/p.fs_in;

% Interpolation time points in samples
if p.tx_angle >= 0
    t_ix    = p.fs_in*((cos(p.tx_angle)*Z + X*sin(p.tx_angle) + sqrt(Z.*Z + (X_ch - X).*(X_ch - X)))./p.c - p.t0);
else
    t_ix    = p.fs_in*((cos(p.tx_angle)*Z + ((p.nch-1)*p.dx - X)*sin(-p.tx_angle) + sqrt(Z.*Z + (X_ch - X).*(X_ch - X)))./p.c - p.t0);
end

dec = p.fs_in/p.fs_out;

switch type
    case 'nearest'
        % Find nearest sample
        t_ix_floor = round(t_ix);
        apod       = apod.*(t_ix_floor > 0 && t_ix_floor < p.nsamples);
        z_ix       = min(max(t_ix_floor,1),p.nsamples);        
        
        % Compute matrix indices, linear. See sub2ind
        ind    = ((z_ix+1) + (ch_ix-1)*p.nsamples);        
        
        % Compute the beamformed signal
        for nn=1:p.nframes
            bf(:,:,nn) = sum(...
                apod.*ch_data(ind).*exp(1i*w_demod*(z_ix - dec*(iz-1)))...
                ,3);
        end

    case 'linear'       
        % Interplation point lies in between y1 and y2. z_ix is the index
        % to y1 and z_ix_frac is mu, the subsample shift
        % Find linear interpolation coefficients
        t_ix_floor = floor(t_ix);        
        apod       = apod.*(t_ix_floor > 0 & t_ix_floor <= p.nsamples);
        z_ix       = min(max(t_ix_floor,1),p.nsamples-1);
        z_ix_frac  = t_ix - z_ix;
            
        % Compute matrix indices, linear. See sub2ind
        ind        = ((z_ix+1) + (ch_ix-1)*p.nsamples);
        
        % Compute the beamformed signal. 
        % x = y1 + (y2 - y1)*mu
        
%         for nn=1:p.nframes
%             bf(:,:,nn) = sum(...
%                 apod.*(ch_data(ind) + (ch_data(ind+1) - ch_data(ind)).*z_ix_frac)...
%                 .*exp(1i*w_demod*(t_ix - dec*(iz-1)))...
%                 ,3);
%         end        
        bf = sum(...
                apod.*(ch_data(ind) + (ch_data(ind+1) - ch_data(ind)).*z_ix_frac)...
                .*exp(1i*w_demod*(t_ix - dec*(iz-1)))...
                ,3);
    case 'cubic' 
        % http://paulbourke.net/miscellaneous/interpolation/
        
        % Interplation point lies in between y1 and y2. z_ix is the index
        % to y1 and z_ix_frac is mu, the subsample shift
        % Find linear interpolation coefficients
        t_ix_floor = floor(t_ix);
        apod       = apod.*(t_ix_floor > 0 && t_ix_floor < p.nsamples-1);
        z_ix       = min(max(t_ix_floor,1),p.nsamples-2);
        z_ix_frac  = t_ix - z_ix;

        % Compute matrix indices, linear. See sub2ind
        ind        = ((z_ix-1)+1 + (ch_ix-1)*p.nsamples);

        % Compute the beamformed signal
        % x = y1 + (-y0 + y2)*mu + (2*y0 - 2*y1 - y3 + y2)*mu^2 + (-y0 + y1
        % - y2 + y3)*mu^3
        for nn=1:p.nframes
            bf(:,:,nn) = sum(...
                apod.*( (...
                (-ch_data(ind-1) + ch_data(ind) - ch_data(ind+1) + ch_data(ind+2)).*z_ix_frac + ...
                (2*ch_data(ind-1) - 2*ch_data(ind) + ch_data(ind+1) - 2*ch_data(ind+2)).*z_ix_frac +...
                (-ch_data(ind-1) + 2*ch_data(ind+1))).*z_ix_frac +...
                ch_data(ind))...
                .*exp(1i*w_demod*(t_ix - dec*(iz-1)))...
                ,3);
        end   
end

% Define the interpolation grids
X = X(:,:,1);
Z = Z(:,:,1);


function apod = gen_apodization(apod,apod_type,varargin)
    
    switch apod_type
        case 'rect'
            % Do nothing
            apod = double(apod >= 1);
        case 'hamming'
            apod = (apod >= 1).*(0.54 + 0.46*cos(pi./apod));             
        case 'tukey'           
            apod = double(apod >= 1);
            for iz=1:size(apod,1)
                for ix=1:size(apod,2)
                    ii = find(apod(iz,ix,:) >= 1);
                    apod(iz,ix,ii) = tukeywin(length(ii),varargin{1});
                end
            end            
        case 'cheb'
            apod = double(apod);
            for iz=1:size(apod,1)
                for ix=1:size(apod,2)
                    ii = find(apod(iz,ix,:));
                    apod(iz,ix,ii) = chebwin(length(ii),varargin{1});
                end
            end            
    end


