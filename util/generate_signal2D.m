function sig = generate_signal2D(u,N,depths,z_mm,win_h,normalize,hr)

Nfft = size(N,1);
Nbeams = size(N,2);

Nt = size(u,1);
Naz = size(u,2);

% Weights
w_total = zeros(Nfft,1);

% Output signal
sig = zeros(Nfft,Nbeams);

x_offset = floor((Nfft - Nt)/2);
y_offset = floor((Nbeams - Naz)/2);

% fh = figure(100);
% clf(fh)
alpha = 0.6;
%cmap = colormap(lines(length(depths)));
for kk=2:length(depths)    
    z_kk = depths(kk);
    
    u_kk = u(:,:,kk);        
    if normalize
        u_kk = u_kk./max(u_kk(:));
    end
    
    U_kk = fft2(hr(:,:,kk).*u_kk,Nfft,Nbeams);
    
    SIG_kk = U_kk.*N;
    sig_kk = fftshift(real(ifft2(SIG_kk)));
    sig_kk = circshift(sig_kk,[x_offset y_offset]);
    
    z_ix = find_indx(z_mm,z_kk);
    
    win_ix = z_ix + (-win_h:win_h);
    win_ix = win_ix(win_ix <= Nfft);
    
    win_kk = zeros(Nfft,1);
    win_kk(win_ix) = tukeywin(length(win_ix),alpha);
    
    w_total = w_total + win_kk;
    sig = sig + win_kk(:,ones(Nbeams,1)).*sig_kk;
    
%     figure(100)
%     subplot(211)
%     hold on;
%     plot(z_mm,win_kk);    
%     set(ph,'Color',cmap(kk,:));
%     hold off
%     subplot(212)
%     hold on;
%     ph = plot(z_mm(win_ix),win_kk(win_ix,ones(Nbeams,1)).*sig_kk(win_ix));
%     set(ph,'Color',cmap(kk,:));
%     hold off
    
%     pause;
end
w_total(w_total == 0) = 1;
sig = sig./w_total(:,ones(Nbeams,1));
