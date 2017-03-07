function make_analytical_signal(h)
%MAKE_ANALYTICAL_SIGNAL Summary of this function goes here
%   Detailed explanation goes here

assert(h.format==E.signal_format.RF,'Signal format is not RF. Analytical signal can only be made from RF data.');

h.data = hilbert(h.data);
h.format = E.signal_format.AS;
end

