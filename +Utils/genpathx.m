function p = genpathx(d,varargin)
% GENPATHX generate recursive path while excluding certain patterns.
%
%     p = GENPATHX(d,patterns) returns a path string starting in d, plus,
%     recursively, all the subdirectories of d, including empty
%     subdirectories, but excluding folders that match the patterns in
%     'pattern'. Each pattern can contain wildcards.
%
%     p = GENPATHX(d,patterns,'-regexp') uses regular expressions to match
%     patterns.
%
%     Examples:
%       p = GENPATHX(d,'.svn');
%       p = GENPATHX(d,{'.svn','dummy*'});
%       p = GENPATHX(d,'^\.\w+$','-regexp');
%

if nargin==0,
    p = genpath(fullfile(matlabroot,'toolbox'));
    if length(p) > 1, p(end) = []; end % Remove trailing pathsep
    return
end

% --------------------
if nargin == 1,
    p = genpath(d);
    return
end
opt = [];
x = cellstr(varargin{1}); % make sure it's a cell string.
if nargin < 3 && ~any(strcmpi(opt,{'-r','-regexp'})),
    x = cellfun(@(y)['^',regexptranslate('wildcard',y),'$'],x,'un',0);
end
% >= R2012b has 'strjoin' command
try
    str = [sprintf('%s|',x{1:end-1}),x{end}];
catch
    str = '';
end
% --------------------

% initialise variables
classsep = '@';  % qualifier for overloaded class directories
packagesep = '+';  % qualifier for overloaded package directories
p = '';           % path to be returned

% Generate path based on given root directory
files = dir(d);
if isempty(files)
    return
end

% Add d to the path even if it is empty.
p = [p d pathsep];

% set logical vector for subdirectory entries in d
isdir = logical(cat(1,files.isdir));
%
% Recursively descend through directories which are neither
% private nor "class" directories.
%
dirs = files(isdir); % select only directory entries from the current listing

for i=1:length(dirs)
    dirname = dirs(i).name;
    if    ~strcmp( dirname,'.')          && ...
            ~strcmp( dirname,'..')         && ...
            ~strncmp( dirname,classsep,1) && ...
            ~strncmp( dirname,packagesep,1) && ...
            ~strcmp( dirname,'private') && ...
            isempty(regexp(dirname,str,'once'))
        p = [p genpathx(fullfile(d,dirname),x,opt)]; % recursive calling of this function.
    end
end