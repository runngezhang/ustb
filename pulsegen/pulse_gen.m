function p = pulse_gen(fc,Ts,Nc,amp,form,envform)

if nargin < 4
    form = 'sin';
    envform = 'none';
    amp = 1;
elseif nargin < 5
    form = 'sin';
    envform = 'none';
elseif nargin < 6
    envform = 'none';
end

Tp = Nc/fc;

N = ceil(Tp/Ts) + 1;
t = Ts*(0:(N-1))';

env = envelope(fc,Nc,N,Ts,envform);

switch form
    case 'sin'
        p = amp*env.*sin(2*pi*fc*t);
    case 'cos'
        p = amp*env.*cos(2*pi*fc*t);
    case 'square'
        p = amp*env.*square(2*pi*fc*t);
end 