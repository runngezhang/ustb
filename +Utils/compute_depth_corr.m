function corr = compute_depth_corr(x,win)

N = size(x,1);
M = size(x,2);
corr = ones(N,M);

for kk=(1+win):(N-win)
    
    x_kk = x(kk + (-win:win),:,:);
    
    x_11 = (sum(x_kk(:,:,1).*conj(x_kk(:,:,1)),1));
    x_22 = (sum(x_kk(:,:,2).*conj(x_kk(:,:,2)),1));
    
    x_12 = sum(x_kk(:,:,1).*conj(x_kk(:,:,2)),1);
    
    corr(kk,:) = real(x_12./sqrt(x_11.*x_22));
end