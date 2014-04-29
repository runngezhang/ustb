function mip_seq = make_mip(img_seq,start_frame,end_frame,roi,lambda,init_frame)
% mip_seq = make_mip(img_seq,start_frame,end_frame,roi)
% 
% Function for generation maximum intensity projection (MIP) of an image
% sequence. The maximum intensity projection is:
%           
%           scale = 1 - lambda;
%           mip_img[1] = init_frame;
%           for k=1:Nframes
%               mip_img[k] = max(scale*mip_img[k-1],img[k]);
%           end
%
% Input:
% img_seq       - Image sequence
% start_frame   - (Optional) Indx of first frame. If set to [], defaults to
%                  1.
% end_frame     - (Optional) Indx of last frame. If set to [], defaults to
%                  all frames
% roi           - (Optional) Region of interest, (x0,z0,x1,z1) where x is
%                 lateral dimension and z is depth dimension. If set to
%                 [], the whole region is used.
% lambda        - (Optional) Used to scale the previous MIP image by
%                 (1-lambda) before doing comparison with current image.
%                 Default = 0.
% init_frame    - (Optional) The frame to initialize the mip generation
%                 with. Default img_seq[start_frame]

% Get dimensions of input
dimensions  = size(img_seq);
Nframes     = dimensions(end);
Ndepth      = dimensions(1);
Nwidth      = dimensions(2);

% Check if start and end frame indxs are provided
if nargin < 2 || isempty(start_frame)
    start_frame = 1;
end

if nargin < 3 || isempty(end_frame)
    end_frame = Nframes;
else
    end_frame = min(Nframes,end_frame);
end

if nargin < 4 || isempty(roi)
    z_indx = 1:Ndepth;
    x_indx = 1:Nwidth;        
else
    z_indx = roi(2):min(roi(4),Ndepth);
    x_indx = roi(1):min(roi(3),Nwidth);    
end

if nargin < 5 
    lambda = 0;
end
scale = 1-lambda;

if nargin < 6
    init_frame = img_seq(:,:,start_frame);    
end

% Calculate the correct number of frames
frame_indxs = start_frame:end_frame;
Nframes = length(frame_indxs); 

% Remove singleton dimensions
img_seq = squeeze(img_seq);

% Allocate mip_sequence matrix
mip_seq = zeros(Ndepth,Nwidth,Nframes,class(img_seq));

% Set the first frame
mip_seq(z_indx,x_indx,1) = init_frame(z_indx,x_indx);

% Calculate MIP sequence
for k=2:Nframes    
    mip_seq(z_indx,x_indx,k) = max(scale*double(mip_seq(z_indx,x_indx,k-1)),double(img_seq(z_indx,x_indx,frame_indxs(k))));            
end