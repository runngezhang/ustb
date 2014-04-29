function [s,r] = pdflatex(texfile)
% Compile latex file with pdflatex
% pdflatex myfolder\myfile.tex (extension .tex is optional)

[fpath,filename,ext] = fileparts(texfile);
[s,r] = system(sprintf('pdflatex -output-directory=%s %s',fpath,[filename,ext]));
if s && ~nargout,
    fprintf('%s\n',r);
end