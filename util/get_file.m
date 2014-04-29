function filename = get_file(filepath,extension,selectFile)
% filename = get_file(filepath,extension)
%
% Lists the files with the specified extension in the specified filepath
% and lets you choose one of them. The full filepath is returned. If no
% extension is specified, '.rf' is used. If cancel is clicked the filename
% is empty.
%
% Input:
% filepath      - path to list files from
% extension     - the filetype you want listed
% selectFile    - optional: if filename (without path) given, not empty and 
%                 present in the given directory, it will be selected in 
%                 the list
%
% Output:
% filename      - the full path to the choosen file.

if nargin == 1
    extension = '.rf';
end

if ~exist(filepath,'dir')
    error(['Directory does not exist : ' filepath]);
end
if ~exist('selectFile','var')
    selectFile = '';
end;
checkSelectFile = ~isempty(selectFile);
selectedIdx = [];

contents = dir(filepath);
files = [];
cnt = 1;
for k=1:length(contents)
    if ~contents(k).isdir
        [dummy, name, ext] = fileparts(contents(k).name);
        if strcmp(ext,extension)
            files{cnt} = contents(k).name;
            if checkSelectFile && strcmp(files{cnt},selectFile)
                selectedIdx = cnt;
            end;
            cnt = cnt + 1;
        end
    end
end

if isempty(files)
    error(['No ' extension ' files in ' filepath]);
end

if isempty(selectedIdx)
    [selection ok] = listdlg('ListString',files,'SelectionMode','single');
else
    [selection ok] = listdlg('ListString',files,'SelectionMode','single','InitialValue',selectedIdx);
end

if ~ok
    filename = [];
else
    filename = fullfile(filepath,files{selection});
end