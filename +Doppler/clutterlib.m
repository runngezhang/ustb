function res=clutterlib(ftype,Iwant,N,Nb,Nf,Fperiod,attenatzero,m,w0);
%CLUTTERLIB  -   utility functions for color flow clutter filtering
%res=clutterlib(ftype,Iwant,N,Nb,Nf,Fperiod,attenatzero,ms);
%ftype: 'poly' (Legendre polynomials),
%       'fourier' (bk(n)=exp(i*2*pi*k*n/Fperiod), sorted by increasing abs.
%       frequency
%Iwant: 'base', 'filtermatrix','freqrespons'
%N: packet size, N>2
%Nb: number of removed components for regression filters
% For FIR, Nb is the filter coefficients (row vector)
%Nf: Number of frequency points in freqrespons over -pi..pi, default=1024
%Fperiod: period for Fourier basis functions, default = N
%attenatzero: attenuation at f=0 i dB. Default attenatzero=10000
% m is the lag in autocorrelation function R(m) for freq. response calculation 

% Author Hans Torp  31.08.04
% Included abs(Hm/H0) plot in freqrespons, 03.07.05  HT
% Included zero at f=0 when abs(w0)>0  11.03.09  HT
% Larger, skewed stopband; Nb= 2D vector allowed  11.03.09  HT

if nargin<1, 
    help clutterlib;return;
end;
if nargin<4, Nb=2;end;
if nargin<5, Nf=1024;end;
if nargin<6, Fperiod=N;end;
if nargin<7, attenatzero=10000;end;
if nargin<8, m=0;end;
if nargin<9, w0=0;end;
base=[];
switch lower(ftype),
    case 'poly', 
        base=polybase(N,Nb,w0);
        base=GrammSchmith(base);
    case 'fourier',
        base=fourierbase(N,Fperiod,-mod(Nb-1,2)/2);
        base=GrammSchmith(base);
    case 'fir',
        base=eye(N);
        if length(Nb)==1, Nb=fir1(Nb,0.2,'high');end;
    otherwise disp('Unknown ftype');return;
end;

switch lower(Iwant),
    case 'base', 
        res=base';
    case 'filtermatrix',
        res=regfiltermatrix(base,Nb,attenatzero); 
    case {'freqrespons','fresp'}
        Fm=regfiltermatrix(base,Nb,attenatzero); 
        [H,f]=matfresp(Fm',Nf,m); 
        [H0,f]=matfresp(Fm',Nf,0); 
        if nargout<1, dispfresp(H,H0,f); else res=H;end;
    otherwise disp('What do you really want from CLUTTERLIB ?');
end;

function base=polybase(N,Nb,w0);
t0=linspace(-0.5,0.5,N);
t0=t0';
base=zeros(N,N);
if abs(w0)==0|length(Nb)<2,
    for p=0:N-1,
        base(:,p+1)=t0.^p;
    end;
else
    for p=0:Nb(1)-1,
        base(:,p+1)=t0.^p ;
    end;
    for p=Nb(1):N-1,
        %base(:,p+1)=t0.^(p-Nb(1)) .*exp(i*w0*t0*N);
        base(:,p+1)=t0.^0 .*exp(i*(p-Nb(1)+1)^0.63*w0*t0*N);
    end;
end;
function base=fourierbase(N,Nf,df);
%Fourier basis, N=packetsize, Nf= interval length
% df=frequency offset; use df=-0.5 for odd order
k=0:N-1;
n=0;
for nn=1:N,
    n=[n;nn;-nn];
end;
n=n(1:N)+df;
base=exp(i*2*pi/Nf*n*k);
base=conj(base');

function base=GrammSchmith(base0);
%Gramm-Schmith orthonormalization
base0=base0';
N=size(base0,1);
base=zeros(N,N);
for p=0:N-1,
    b=base0(p+1,:);
	for k=0:p-1, %ortogonalization
			 bk=base(k+1,:);
			 b=b-(b*bk')*bk;
	end;
	b=b/sqrt(b*b');%Normalization
	base(p+1,:)=b;
end;%for p
base=base';

function Fm=regfiltermatrix(base,Nb,attenuationatzero);
% regression filter matrix 
% Fm:  filter matrix
% base: matrix with basis functions
% For FIR-filters,  base is a row vector of filter coefficients
% Nb: #basisfunctions to remove
% w0 center frequency of stopband
if length(Nb)>1, Nb=sum(Nb);end;
N=size(base,1);
b0=ones(N,1)/sqrt(N);
if length(Nb)>1, Fm=firmatr(Nb,N);
else
    Fm=eye(N);
    if Nb>0&Nb<N,
        bp=base(:,1:Nb);
        c=1-10^(-attenuationatzero/20);
        Fm = Fm - c*bp*bp';
    end;
end;


function [H,f]=matfresp(Fm,Nf,m);
% [H,f]= matfresp(Fm,m,Nf);
%calculate frequency response for R(m) for clutterfilter given by the filter matrix Fm
%m=temp.lag in correlation function
%Fm=filter matrix, NB: Input when using matfresp is Hermetian transpose of Fm!!!

% Author     :   Hans Torp 4/2-97
%Fm=Fm';
[M,N]=size(Fm);
f=0:1/Nf:(1-1/Nf); 
w=2*pi*f';
f=f-0.5;
A=fft(Fm,Nf);
Ap=conj(A(:,1:N-m)).*A(:,1+m:N);
H=fftshift(mean(Ap').*exp(i*m*w)');

function [H,f]=matfresp0(Fm,Nf,m);
H=mean(abs(fft(Fm,Nf)).^2,2);
f=0:1/Nf:(1-1/Nf); 

function dispfresp(H,H0,f);
%clf;
plot(f,10*log10(abs(H)));grid;axis([-0.5 0.5 -80 0]);
ylabel('Amplitude [dB]');
if ~isreal(H)
    subplot(3,1,1);	plot(f,10*log10(abs(H0)));grid;axis([-0.5 0.5 -80 0]);
    ylabel('H0 Amplitude [dB]');
    subplot(3,1,2);	plot(f,abs(H./(H0+1e-30)));grid;ylabel('abs(Hm/H0)');
    subplot(3,1,3);	plot(f,angle(H)/pi/2);grid;ylabel('Bias frequency/fs');
end;
xlabel('Frequency/fs');
function fmatr=firmatr(h,psize);
%filtermatrix for FIR filter
order=length(h)-1;
h=fliplr(h);
fmatr=zeros(psize-order,psize);
for n=1:psize-order
	fmatr(n,n:n+order)=h;
end;
