function savefig(fh,fig_path,fname,do_save)

if nargin < 4
    do_save = 0;
end

if do_save    
    %saveas(fh,fullfile(fig_path,fname),'png');
    saveas(fh,fullfile(fig_path,fname),'fig');
    %saveas(fh,fullfile(fig_path,[fname,'.eps']),'psc2');    
    export_fig(fullfile(fig_path,fname),'-eps','-png','-pdf','-nocrop',fh);
end