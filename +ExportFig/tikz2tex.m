function tikz2tex(filename,preview)
% Wrap tikz/pgf in latex-environment
% tikz2tex(filename,preview)
%
[fpath,fname,ext] = fileparts(filename);
texfilename = fullfile(fpath,[fname,'.tex']);
pdffilename = fullfile(fpath,[fname,'.pdf']);

fid = fopen(texfilename,'w+');
% write tex-file
fprintf(fid,'%s\n','\documentclass{figview}');
fprintf(fid,'%s\n','\usepackage{tikz,siunitx}');
fprintf(fid,'%s\n','\usepackage[active,tightpage]{preview}');
fprintf(fid,'%s\n','\PreviewEnvironment{tikzpicture}');
fprintf(fid,'%s\n','\setlength\PreviewBorder{1pt}%');
fprintf(fid,'%s\n','\begin{document}');
fprintf(fid,'\\input{%s}\n',strrep([fname,ext],'\','/'));
fprintf(fid,'%s\n','\end{document}');
fclose(fid);

% compile latex and show pdf file
pdflatex(texfilename);
if nargin > 1 && preview,
    system(sprintf('SumatraPDF %s &',pdffilename));
end