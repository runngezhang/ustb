function Fm = polyregMatrix(N,P)
% polyregMatrix creates filter matrix for polynomial regression filter
% Fm=polyregMatrix(N,P) Generate NxN filter matrix Fm for polynomial regression filter
% N = packet size (number of temporal samples)  N >= 2
% P = Polynomial order .  OBS:  -1<= P < N-1.  P=-1 gives Fm=identity matrix

% Author :   Hans Torp 15.10.01

Fm= zeros(N,N);
if P> (N-2), disp('Filter order too large');return; end;
if P<0, Fm = eye(N,N);return;end;%identity matrix

%form the polynomial basis vectors by Gramm-Schmith
t=1:N;
base=zeros(P+1,N);
for p=0:P,
	b=t.^p;%polynom of order p
	for k=0:p-1, %ortogonalization
		bk=base(k+1,:);
		b=b-(b*bk')*bk;
	end; 
	b=b/sqrt(b*b');%Normalization
	base(p+1,:)=b;
	Fm=Fm+b'*b;
end;%for p
Fm = eye(N,N) - Fm;
