function [result, modifiedFiles] = checkVCstate(functionList, fileList, zipFilename)
    % [result, modifiedFiles] = checkVCstate(functionList, fileList, zipFilename)
    %
    % checks the version control state: 
    %    - under version control?
    %    - locally modified?
    %    - last committed revision
    %    - updated revision
    % of the given files (cell arrays: fileList) and functions (it will be
    % determined before checking the status which file is called by a
    % certain function name)
    %
    % The results are returned in an array of a structure: 
    % result.filename
    %       .underVC
    %       .isLocallyModified
    %       .lastCommitRev
    %       .updatedRev
    % modifiedFiles contains a list with the filenames of files which are
    %       modified or not under version control.
    %       
    % If zipFilename is given all functions with local modifications or 
    % which are not under version control will be zipped together in this 
    % file.
    % Author: Jochen Deibele, 2010
    
    vcpath = 'SubWCRev.exe';
    % check if tool is present
    output = evalc(sprintf('!%s',vcpath));
    lines = regexp(output,'\n','split');
    if isempty(regexp(lines{1},'SubWCRev', 'once'))
        error('Subversion version checker SubWCRev not found');
    end;
    
    % create empty result structure
    result = struct('filename',{},'isLocallyModified',{},'lastCommitRev',{},...
                    'updatedRev',{},'underVC',{});
    
    % we start with no modified files
    modifiedFiles = {};
    
    % find file names for functions
    for i = functionList
        fcn = i{1}; 
        
        % get full filename
        filename = evalin('caller',sprintf('which(''%s'')',fcn));
        if isempty(filename)
            warning('function %s not found',fcn);
        else
            fileList{end+1}= filename;
        end;
    end
    
    % check status of files
    fprintf('\ncheckVCstate -- checking:\n');
    for i = fileList
        filename = i{1};
        
        if ~exist(filename,'file')
            warning('file %s not found',filename);
        else
            cur.filename = filename;
            msg = sprintf('   * %s ...',filename);
            fprintf('   * %s ...',filename);
            % has local modifications?
            [hasLocalMod, res] = system(sprintf('%s %s -n',vcpath,filename));
            cur.isLocallyModified = hasLocalMod == 7;

            % get revision        
            output = evalc(sprintf('!%s %s',vcpath,filename));
            % parse output for revision
            lines = regexp(output,'\n','split');

            cur.lastCommitRev = 0;
            cur.updatedRev = 0;
            cur.underVC = false;
            if numel(lines)>2
                lastCommitRev = regexp(lines{2},'Last committed at revision (\d+)','tokens');
                lastCommitRev = lastCommitRev{1}; % unwrap cell arrays
                cur.lastCommitRev = str2double(lastCommitRev{1});
                updatedRev = regexp(lines{3},'Updated to revision (\d+)','tokens');
                updatedRev = updatedRev{1};
                cur.updatedRev = str2double(updatedRev{1});
                % state if file is under version control
                cur.underVC = true;
                if cur.updatedRev == 0 && cur.lastCommitRev == 0
                    cur.underVC = false;
                end;
            end;
            if cur.isLocallyModified || ~cur.underVC
                modifiedFiles{end+1} = cur.filename;
                if cur.isLocallyModified
                    fprintf('   DIRTY\n');
                else
                    fprintf('   NOT VERSIONED\n');
                end;
            else
                astr = '\b';
                msg = repmat(astr,1,numel(msg)); % back to start
                fprintf(msg);
                %fprintf('   clean\n');                
            end;            
            result(end+1) = cur;
        end;        
    end;
    fprintf('\n\n');
    
    % if zipfilename given, pack files
    if exist('zipFilename','var') && ~isempty(zipFilename)
        fprintf('the following files were modified or not under version control they''ll be packed:\n');
        for i = modifiedFiles
            fprintf('* %s\n',i{1});
        end;
        zip(zipFilename,modifiedFiles);
    end;
end