function ustb_data = ustb_data_path()
%USTB_PATH Returns the path (root) of data for USTB.
% Returns the value of environment variable 'USTB_DATA' if it exists,
% otherwise it returns path to `data` folder in USTB
%
%   Usage: out = ustb_data_path()
%
    ustb_data = getenv('USTB_DATA'); 
    
    if isempty(ustb_data),
        [ustb_path, ~, ~] = fileparts(which('ustb_data_path'));
        ustb_data = [ustb_path, filesep, 'data'];
    end
    
    if ~exist(ustb_data, 'dir'),
        warning('USTB:DataPathNotExist', ['"', ustb_data, '"', ' does not exist!']);
    end
end
