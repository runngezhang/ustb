function fc = freq_estimation(x,fs,win_h,step,method,Nfft)
% fc = freq_estimation(x,fs,win_h,step,method)
% 
% Function for estimating center frequency in signal traces
% 
% Input:
% x         - Signal, IQ or RF
% fs        - Sampling frequency
% win_h     - One sided window length to use for estimation
% step      - Depth increment, number of samples between each estimated point.
% The values in between are interpolated
% method    - Estimation method. For IQ, use 'ar1' and for RF use 'ar2'
%
% Output:
% fc        - Center frequency estimates
%
% For 'ar1' the lag one autocorrelation estimator (Doppler estimator by
% Kasai[1]) is used. For 'ar2' the AR2 estimator Girault[2] is used.
%
% [1]  C. Kasai, K. Namekawa, A. Koyano, and R. Omoto, ?Real-time
% two-dimensional blood flow imaging using an autocorrelation technique,?
% IEEE Trans. Sonics Ultrason. SU-32(3), 458?463 (1985)  
%
% [2] J-M. Girault, F. Ossant, A., D. Kouame, and F. Patat, "Time-Varying
% Autoregressive Spectral Estimation for Ultrasound Attenuation in Tissue
% Characterization" IEEE Trans on Ultrasonics, Ferroelectrics, and
% Frequency control, Vol. 45 (3), 650-659 (1998)

Ns = size(x,1);
Nb = size(x,2);
Nf = size(x,3);

%bw = zeros(Ns,Nb,Nf);
theta_n = [1 1]';
P_n = eye(2);

if strcmp(method,'fft')
    %fc = zeros(Ns,Nb*Nf);
    fc = zeros(Ns,1);
    
    for c_ix=[1:step:Ns Ns]
        
        ixs = c_ix + (-win_h:win_h);
        ixs = max(1,min(ixs,Ns));

        N = length(ixs);
        if nargin < 6
            Nfft = N;
        end

        win = hamming(N);                
        
        xx = x(ixs,:);
        win = win(:,ones(1,size(xx,2)));
        
        Xf = abs(fft(win.*xx,Nfft,1));        
        Xf = mean(Xf,2);
        
        f = fs*(0:(Nfft-1))'/Nfft;

        Nfft2 = floor(Nfft/2);
        df = f(2) - f(1);

        P = sum(Xf(1:Nfft2,:)*df,1);
        fc(c_ix,:) = (1./P).*sum(Xf(1:Nfft2,:).*f(1:Nfft2,ones(1,size(Xf,2)))*df);                               
    end
    
    d_ixs = 1:Ns;
    e_ixs = 1:step:Ns;
    i_ixs = setdiff(d_ixs, e_ixs);
    fc(i_ixs,:) = interp1(e_ixs',fc(e_ixs,:),i_ixs','linear','extrap');    
    %fc = reshape(fc,[Ns Nb Nf]);
else
    fc = zeros(Ns,Nb,Nf);
    
    for ii=1:Nf
        for kk=1:Nb                        
            for c_ix=[1:step:Ns Ns]

                ixs = c_ix + (-win_h:win_h);

                switch method
                    case 'ar1'

                        ixs = max(1,min(ixs,Ns-1));

                        if isreal(x)
                            error('For ar1 method, IQ signal must be used')
                        end

                        R1 = mean(x(ixs + 1,kk,ii).*conj(x(ixs,kk,ii)));
                        fc(c_ix,kk,ii) = fs*angle(R1)/2/pi;

    %                     if nargout > 1                        
    %                         R0 = mean(x(ixs,kk,ii).*conj(x(ixs,kk,ii)));
    %                         bw(c_ix,kk,ii) = fs*sqrt(2)*sqrt(1 - min(abs(R1/R0),1))/2/pi;
    %                     end

                    case 'ar2'
                        ixs = max(1,min(ixs,Ns-2));

                        if ~isreal(x)
                            error('For ar2 method, RF signal must be used')
                        end                          
                        N = length(ixs);

                        R0 = (1/N)*(sum(x(ixs,kk,ii).*x(ixs,kk,ii)));
                        R1 = (1/N)*(sum(x(ixs+1,kk,ii).*x(ixs,kk,ii)));
                        R2 = (1/N)*(sum(x(ixs+2,kk,ii).*x(ixs,kk,ii)));

                        R_mtrx = [R0 R1; R1 R0];
                        r_mtrx = [-R1; -R2];

                        a = inv(R_mtrx)*r_mtrx;
                        a1 = a(1);
                        a2 = a(2);

                        fc(c_ix,kk,ii) = real(fs/(2*pi)*acos( (-a1/4) *( 1 + 1/a2 ) ));
                    case 'fft'
                        ixs = max(1,min(ixs,Ns));

                        N = length(ixs);
                        if nargin < 6
                            Nfft = N;
                        end

                        win = hamming(N);
                        Xf = abs(fft(win.*x(ixs,kk,ii),Nfft));


                        f = fs*(0:(Nfft-1))'/Nfft;

                        Nfft2 = floor(Nfft/2);
                        df = f(2) - f(1);

                        P = sum(Xf(1:Nfft2)*df);

                        fc(c_ix,kk,ii) = (1/P)*sum(Xf(1:Nfft2).*f(1:Nfft2)*df);

                    case 'ar2seq'                
                        lambda = 0.98;
                        if ~isreal(x)
                            error('For ar2 method, RF signal must be used')
                        end                          
                        theta_n1 = theta_n;
                        P_n1 = P_n;

                        phi_n = -x(c_ix + [-1 -2],kk);
                        eps_n = x(c_ix) - phi_n'*theta_n1;
                        P_n = (P_n1/lambda)*(1 - (P_n1*(phi_n*phi_n')*P_n1')/(lambda + phi_n'*P_n1*phi_n));
                        theta_n = theta_n1 + P_n*phi_n*eps_n;

                        a1 = theta_n(1);
                        a2 = theta_n(2);
                        fc(c_ix,kk) = real(fs/(2*pi)*acos( (-a1/4) *( 1 + 1/a2 ) ));
                        c_ix = c_ix + 1;
                end           
            end

            d_ixs = 1:Ns;
            e_ixs = 1:step:Ns;
            i_ixs = setdiff(d_ixs, e_ixs);
            fc(i_ixs,kk,ii) = interp1(e_ixs,fc(e_ixs,kk,ii),i_ixs,'linear','extrap');    
     %       bw(i_ixs,kk,ii) = interp1(e_ixs,bw(e_ixs,kk,ii),i_ixs,'linear','extrap');    
        end
    end
end