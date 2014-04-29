function [n_bp,h,N_bp,H,F] = bandpass_noise(Nfft,Nseqs,fs,fco,Ntaps,silent,win_type,domain)
% bandpass_noise(Nfft,fs,fco,Ntaps,silent,win_type)
%
% Function for generating bandpass filtered noise. The filtering is done in
% the FFT domain. The frequency response of the filter is computed using
% freqz. The filtering can also be done in the time domain by specifying
% this.
%
% Input : 
% Nfft      - Length of noise sequence
% Nseqs     - Number of sequences
% fs        - Sampling frequency
% fco       - Cutoff frequencies of BP filter
% Ntaps     - Number of taps of BP filter
% silent    - (Optional) If 1 the spectrum of the sequence and the filter
% is plotted. Default = 0
% win_type  - (Optional) Window for the BP filter. hamming, gauss, tukey,
% ones. Default = 'hamming'. See fir1
% domain    - (Optional) In what domain the filtering should be done, 'fft'
% or 'time'. Default = 'fft';
%
% Output
% n_bp      - The filtered noise sequence
% h         - BP filter
% N_bp      - Sepctrum of filtered noise
% H         - Frequency response of h, generated using freqz
% F         - Frequency axes

if nargin < 1
	Nfft = 10e3;
end
	
if nargin < 2
	Nseqs = 1;
end

if nargin < 3
	fs = 50e6;
end

if nargin < 4
	fco = [0.75 18]*1e6;
end

if nargin < 5
	Ntaps = 100;
end

if nargin < 6
    silent = 0;
end

if nargin < 7
    win_type = 'hamming';
end

if nargin < 8
    domain = 'fft';
end

Ts = 1/fs;	
switch win_type
    case 'hamming'
        filter_win = hamming(Ntaps+1);
    case 'gauss'        
        filter_win = gausswin(Ntaps+1);
    case 'tukey'        
        filter_win = tukeywin(Ntaps+1);
    case 'ones'        
        filter_win = ones(Ntaps+1,1);
end

h    = fir1(Ntaps,2*fco*Ts,filter_win);
[H,F] = freqz(h,1,Nfft,'whole',fs);

[db_lim] = find_db_limits(H(1:Nfft/2),-6,F(1:Nfft/2));

n = randn(Nfft,Nseqs);
N = fft(n,[],1);

switch domain
    case 'fft'
        N_bp = H(:,ones(Nseqs,1)).*N;
        n_bp = real(ifft(N_bp,[],1));
    case 'time'
        n_bp = filter(h,1,n,[],1);
        N_bp = fft(n_bp,[],1);
end

if ~silent
    figure;
    subplot(221)
    ph = plot(1e-6*F,20*log10(abs(H/max(H)))...
        ,1e-6*F(db_lim(1)*[1 1]),[0 -80],'k:'...
        ,1e-6*F(db_lim(2)*[1 1]),[0 -80],'k:');
    set(gca,'FontSize',16)
    set(ph,'LineWidth',2);
    title('Filter response')
    xlabel('[MHz]')
    ylabel('[dB]')
    xlim(1e-6*[0 fs/2])
    ylim([-80 0])
    grid on;

    subplot(222)
    ph = plot(1e-6*F,20*log10(abs(N/max(N(:))))...
        ,1e-6*F(db_lim(1)*[1 1]),[0 -80],'k:'...
        ,1e-6*F(db_lim(2)*[1 1]),[0 -80],'k:');
    set(gca,'FontSize',16)
    set(ph,'LineWidth',2);
    title('Noise spectrum')
    xlabel('[MHz]')
    ylabel('[dB]')
    xlim(1e-6*[0 fs/2])
    ylim([-80 0])
    grid on;

    subplot(212)
    ph = plot(1e-6*F,20*log10(abs(N_bp/max(N_bp(:))))...        
        ,1e-6*F(db_lim(1)*[1 1]),[0 -80],'k:'...
        ,1e-6*F(db_lim(2)*[1 1]),[0 -80],'k:');
    set(gca,'FontSize',16)
    set(ph,'LineWidth',2);
    title('Filtered noise')
    xlabel('[MHz]')
    ylabel('[dB]')
    xlim(1e-6*[0 fs/2])
    ylim([-80 0])
    grid on;

%     figure(2)
%     subplot(211)
%     ph = plot(t_ms,n,t_ms,n_bp);
%     set(gca,'FontSize',16)
%     set(ph,'LineWidth',2);
%     title('Time signal')
%     xlabel('[\mus]')
%     ylabel('[.]')
%     grid on;
% 
%     c = 1.540;
%     z_mm = c*t_ms/2;
%     subplot(212)
%     ph = plot(z_mm,n,z_mm,n_bp);
%     set(gca,'FontSize',16)
%     set(ph,'LineWidth',2);
%     title('Time signal')
%     xlabel('[mm]')
%     ylabel('[.]')
%     grid on;
end