function make_analytical_signal(h)
%MAKE_ANALYTICAL_SIGNAL Creates the analytical signal by taking the Hilbert
%   transform
%
%   authors : Ole Marius Hoel Rindal <olemarius@olemarius.net>
%             Stine Myren Hverven <stinemhv@ifi.uio.no>
%   $Date   : 2017/03/06 
%             2017/03/09 : Fixed bug for matrices lager than 2D. However
%                          this made the code very slow.

assert(h.format==E.signal_format.RF,'Signal format is not RF. Analytical signal can only be made from RF data.');

wb = waitbar(size(h.data,4),'Taking Hilbert Transform')
for f=1:size(h.data,4)
    waitbar(f/size(h.data,4))
    for m=1:size(h.data,3)
        h.data(:,:,m,f) = hilbert(h.data(:,:,m,f));
    end
end
close(wb)
h.format = E.signal_format.AS;
end

