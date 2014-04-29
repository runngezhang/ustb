% File for clearing all but a group of specified variables. Specify the
% variables you want to keep in your file like this:
% 
% vars_to_keep = {'var1','var2',...};
% clear_all_but;
%
% It is important that the variable is called vars_to_keep. If no
% vars_to_keep is specified only the variable h is kept.

variables = whos;

for k=1:length(variables)
    if strcmp(variables(k).name,'variables') || strcmp(variables(k).name,'vars_to_keep')
        continue;
    end
    
    keep_var = 0;
    
    if exist('vars_to_keep','var')        
        for j=1:length(vars_to_keep)
            if strcmp(vars_to_keep{j},variables(k).name)
                keep_var = 1;
                break;
            end
        end
    elseif strcmp(variables(k).name,'h')
        keep_var = 1;
    end
    
    if ~keep_var
        clear(variables(k).name);
    end       
end

clear('variables','k','vars_to_keep','j','keep_var');
