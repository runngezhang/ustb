% load_into_ws(filename,ws)
%
% Loads the specified .mat file and puts the variables into workspace. If you want to specify a workspace other than 'base', use the second parameter.
%
% Input:
% filename - Name of the .mat file
% ws       - (optional) The workspace where the variables should be put.
function load_into_ws(filename,ws)
h = load(filename);
var_names = fieldnames(h);

if nargin == 1
    ws = 'base';
end

for k=1:length(var_names)        
    assignin(ws,var_names{k},h.(var_names{k}));
end