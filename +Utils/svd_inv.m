function [A_inv Nsv sv U V] = svd_inv(A,th)

[U S V] = svd(A);

sv = diag(S);
sv_rel = cumsum(sv)/sum(sv);
Nsv = find(sv_rel > th,1,'first');

A_inv = V(:,1:Nsv)*diag(1./sv(1:Nsv))*U(:,1:Nsv)';