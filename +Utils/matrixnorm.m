function Anorm = matrixnorm(A)
%Function for normalizing a matrix columnvise
%
%Written by Svein-Erik M?s?y

[col,lines]=size(A);

Amax = ones(col,1)*max(A);

Anorm = A./Amax;
