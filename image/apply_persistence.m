function img_persist = apply_persistence(img_seq,win_length,win_type,init_frame)
% apply_persistence(img_seq,win_length,win_type)
%
% Function for applying persistence. Uses a sliding window over the
% sequence of frames.
%
% Input:
% 
%
if mod(win_length,2) == 0
    warning('Window length must be odd. Increasing window length with 1');
    win_length = win_length + 1;
end

if nargin == 2
    win_type = 'ones';    
end

win_h = (win_length-1)/2;
switch win_type
    case 'ones'
        win = ones(1,win_length);
    case 'gauss'
        win = gausswin(win_length);
    case 'hamming'
        win = hamming(win_length);
    otherwise
        error('Win type must be : ones, gauss or hamming');
end

img_persist = zeros(size(img_seq),class(img_seq));
buffer = zeros(size(img_seq,1),size(img_seq,2));
scale = 0;

if nargin == 4    
    if size(img_seq,1) ~= size(init_frame,1) || size(img_seq,2) ~= size(init_frame,2)
        error('The intial frame must have same size as the input images')
    end
    scale = win_h;
    init_frame = double(init_frame);
else
    init_frame = zeros(size(buffer));
    scale = 0;
end

for kk=1:size(img_seq,3)        
    for ii=(-win_h:0);%win_h)        
        if (kk+ii >= 1) && (kk+ii <= size(img_seq,3))
            buffer = buffer + win(ii+win_h+1)*double(img_seq(:,:,kk+ii));
            scale = scale + win(ii+win_h+1);
        end        
    end
    
    if kk < win_h
        img_persist(:,:,kk) = ((win_h-kk+1)*init_frame + buffer)/scale;
        %scale = (win_h - kk);
        scale=0;
    else
        img_persist(:,:,kk) = buffer/scale;
        scale = 0;
    end
    
    buffer(:) = 0;    
end