function env_enh = image_enhance(iq,LoG_dim,LoG_sigma,do_laplace,G_dim,G_sigma)
% env_enhanced = image_enhance(iq,LoG_dim,LoG_sigma)
%
% Function for image enhancement. Takes the IQ signal and upsamples in the
% lateral direction and averages neighbour lines. Then the envelope is 
% filtered with the laplacian edge enhancement filter and smoothed with a
% gaussian mask. Then then does a Laplacian of Gaussian filtering of the
% envelope as well. The final enhanced envelope is formed by doing and
% average of the 2 filtered envelopes and the original.
%
% Links:
% http://homepages.inf.ed.ac.uk/rbf/HIPR2/log.htm
% http://www.amazon.com/Digital-Image-Processing-Rafael-Gonzalez/dp/020118
% 0758
%
% Input:
% iq            - IQ signal
% LoG_dim       - Dimension of Laplacian of Gaussian mask, range 5-13
% LoG_sigma     - Std deviation of Laplacian of Gaussian mask, range 0.5-2
% do_laplace    - true if laplace edge detection should be done
% G_dim         - Dimension of Gaussian mask, range 5-13 (optional when
% do_laplace = true)
% G_sigma       - Std deviation of Gaussian mask, range 0.5-2 (optional when
% do_laplace = true)
%
% Output:
% env_enh       - Enhanced envelope
%

env = int16(abs(iq));
LoG_mask = fspecial('log',[LoG_dim LoG_dim],LoG_sigma);

env_LoG = abs(imfilter(env,LoG_mask,0,'corr','same'));

if do_laplace    
    
    if nargin == 4
        G_dim = 9;
        G_sigma = 0.8;
    end
    
    G_mask = fspecial('gaussian',G_dim,G_sigma);
    L_mask = fspecial('laplacian',0.3);

    env_L = abs(imfilter(env,L_mask,0,'corr','same'));
    env_LG = imfilter(env_L,G_mask,0,'corr','same');    
    env_enh = uint16(0.33*(env + env_LG + env_LoG));
else
    env_enh = uint16(0.5*(env + env_LoG));
end


