classdef scriptParams < handle
    properties
    
    end
    methods(Static)
        function result = readRawFile(filename,parentParams)
            if ~exist('parentParams','var')
                parentParams = [];
            end;
            fid=fopen(filename);
            content = '';
            while 1
                tline = fgetl(fid);
                if ~ischar(tline), break, end
                content = sprintf('%s%s\n',content,tline);
            end
            fclose(fid);
            parent = parentParams;   %#ok<NASGU> (will be used in eval statement)
            eval(content);
            result = s;
        end
        function msg = printFilename(prefix,filename)
            msg = regexprep(filename, '\\', '/'); % replace backslashes
            msg =['<a href="error:' filename ',3,1">' prefix '</a>' msg]; 
            disp(msg);
            %fprintf('%s\n',msg);
        end
        
        function msg = printFileLinks(params)
            msg = '';
            for i = 1:length(params)
                if mod(i,2)==0
                    msg = [msg '  <a href="error:' params{i} ',3,1">' params{i-1} '</a>  ']; 
                end;
            end;
            disp(msg);
        end
        
        function [params shadowList shadowFieldKeyList] = mergeStructs(child, parent, parentName)
            % merges to structures and keeps an updates shadowlist where
            % the parentName is used. params is the resulting struct.
            shadowList = {};
            shadowFieldKeyList = {};
            params = parent;
            for fieldCell = fieldnames(parent).'
                field = fieldCell{1};
                for keyCell = fieldnames(parent.(field)).'
                    key = keyCell{1};
                    % check if parameter is overwritten locally
                    if isfield(child,field) && isfield(child.(field),key)
%                         msg = sprintf('"%s" is shadowed by %s parameter set\n',key,parentName);
%                         printShadow(msg);
                        shadowList{end+1}=sprintf('%s.%s is shadowed by %s parameter set',field,key,parentName);
                        shadowFieldKeyList{end+1} = sprintf('%s.%s',field,key);
                        params.(field).(key) = child.(field).(key);
                    end;                   
                end
            end
        end
        function printParams(params, shadowList, shadowFieldKeyList)
            for fieldCell = fieldnames(params).'
                field = fieldCell{1};
                fprintf('%s\n',field);
                disp(params.(field))                
                for keyCell = fieldnames(params.(field)).'
                    key = keyCell{1};
                    [ismem, idx] = ismember(sprintf('%s.%s',field,key),shadowFieldKeyList);
                    if ismem
                        fprintf('%s\n',shadowList{idx+3});
                    end;                   
                end
            end
        end
        function [result, globPrmFile, dirPrmFile, filePrmFile,all_params] = loadParameters(setup, directory, filename, quiet)
            if ~exist('quiet','var')
                quiet = 'verbose';
            end;
            if false %exist('cprintf','file')
                printShadow =@(x)cprintf([0.9 0.2 0],sprintf('%s',x));
            else
                printShadow =@(x)fprintf('%s',x);
            end;
            shadowed = false;
                
            %save current path
            currPath = cd;
            
            mfilepath = fileparts(mfilename('fullpath'));
            
            if strcmp(directory,mfilepath)
                result = [];
                error('cannot run utility in the same directory as scriptParams.m');
            end
            shadowList = {};
            if ischar(setup)
                % change directory
                cd([mfilepath,'/scriptParamsSetups'])
                
                setupfilename = sprintf('setup_%s',setup);
                % run function in actual directory
                s_default = feval(setupfilename);
                globalSetupFilename = which(setupfilename);
                globPrmFile = globalSetupFilename;
            else
                s_default = setup; 
                globalSetupFilename = '';
                globPrmFile = '';
            end;
            params = s_default;
            params.basic.globPrmFile = globPrmFile;
            shadowList{end+1} = ...
                scriptParams.printFilename('Glob prm file:', globalSetupFilename);
            
            if isempty(directory)
                directory = pwd;
            end;
            
            % now check for a setupscript in the directory given
            dirParamFilename = sprintf('%s/%s_dir.m',directory,params.basic.setupname);
            shadowList{end+1} = ...
                scriptParams.printFilename('Dir prm file :',dirParamFilename);
            
            % change directory
            cd(directory)
               
            if ~exist(dirParamFilename, 'file')
                disp('creating empty directory parameter file');
                copyfile(sprintf('%s/%s_empty.m',[mfilepath,'/scriptParamsSetups'],params.basic.setupname),dirParamFilename);
            end;
            s_dir = scriptParams.readRawFile(dirParamFilename,params);
            dirPrmFile = dirParamFilename;
            
            % now check for parameter file for file given
            [dummy,f,e] = fileparts(filename);
            if numel(e) > 1 % remove dot from extension
                e = e(2:end);
            end;
            fileParamFilename = sprintf('%s_%s_%s.m',f,params.basic.setupname,e);
            
            fileParamFullFilename = sprintf('%s/%s',directory,fileParamFilename);
            shadowList{end+1} = ...
                  scriptParams.printFilename('File prm file:',fileParamFullFilename);            
            
            shadowFieldKeyList = {};
            % merge params from directory file
            [params shadowListTmp shadowFieldKeyListTmp] = scriptParams.mergeStructs(s_dir, params, 'directory');
            shadowList = [shadowList shadowListTmp];
            shadowFieldKeyList = [shadowFieldKeyList shadowFieldKeyListTmp];
            params.basic.dirPrmFile = dirPrmFile;
            
            

            % change directory
            cd(directory)
            
            if ~exist(fileParamFullFilename, 'file')
                fprintf('creating empty file parameter file\n');
                copyfile(sprintf('%s/%s_empty.m',[mfilepath,'/scriptParamsSetups'],params.basic.setupname),fileParamFullFilename);
            end
            s_file = scriptParams.readRawFile(fileParamFilename,params);
            filePrmFile = fileParamFilename; 
            
            % merge params from file file
            [params shadowListTmp shadowFieldKeyListTmp] = scriptParams.mergeStructs(s_file, params, 'file');
            shadowList = [shadowList shadowListTmp];
            shadowFieldKeyList = [shadowFieldKeyList shadowFieldKeyListTmp];
            params.basic.filePrmFile = filePrmFile;
            
            % print parameters on command line
            if strcmp(quiet,'verbose')
                scriptParams.printParams(params,shadowList, shadowFieldKeyList);
            end;
            result = params;
            % return to 'current' directory
            cd(currPath)
            scriptParams.printFileLinks({'File prm file',fileParamFullFilename, ...
                                         'Dir prm file',dirParamFilename,...
                                         'Glob prm file', globalSetupFilename});            
            if numel(shadowList)>3
                shadowed = true;
            end;
            if shadowed
                for i = 4 : numel(shadowList)
                    printShadow(sprintf('%s\n',shadowList{i}));
                end
            end
            all_params.glob = s_default;
            all_params.dir = s_dir;
            all_params.file = s_file;
        end        
    end
end