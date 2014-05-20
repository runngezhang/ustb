function params = read_paramfile(fname)

fid = fopen(fname);
hcl = onCleanup(@()fclose(fid));

params = struct();

while ~feof(fid)
    line0 = fgetl(fid);
    line = strtrim(regexprep(line0,'#.*','')); % remove trailing comments
    line = regexprep(line,'\(.*?\)',''); % remove stuff between parentheses
    if isempty(line),
        continue
    end
    c = regexp(line,'=','split');
    if numel(c) < 2,
        fprintf('Could not parse line: "%s". Skipping it.\n',line0);
        continue
    end
    paramname = regexprep(strtrim(c{1}),{'[\[\]]','[-/\s]'},{'','_'});
    paramval = str2num(c{end}); %#ok<ST2NM>
    if isempty(paramval), % if it's not numerical, then keep the string.
        paramval = strtrim(c{end});
    end
    params.(paramname) = paramval;
end
