function download(filename, url, local_path)
%DOWNLOAD Downloads a file from an URL to a local path (if needed)
%
% Note: Uses `websave`, which available since Matlab 2014. 
%      For older versions, it requires 'curl' installed.
%      See `https://curl.haxx.se/download.html` for Windows

    if not(exist(local_path, 'dir')),
        mkdir(local_path)
    end
    
    % Make sure that local path ends with file separator
    if ~isempty(local_path) && local_path(end) ~= filesep,
        local_path = [local_path, filesep];
    end
    
    if url(end) ~= '/',
        url = [url, '/'];
    end

    fullpath = [local_path, filename];
    fullurl  = [url,        filename];
    
    if not(exist(fullpath,  'file')),
        disp(['Downloading ', filename, ' from ', url, ' to ', local_path]);
        disp('This may take a while.');

        if exist('websave')
            outfilename = websave(fullpath, fullurl);
            
        elseif find_curl(),
            success = download_with_curl(fullurl, local_path, filename);
            if not(success)
                error(['Could not download ', fullurl, ' to ', fullpath]);
            end
        else
            error('Cannot download data. Install "curl" and try again.');
            
        end
    end
    
end

function [found, curl_cmd] = find_curl()
%Find the executable "CURL" which can be used to download data
%
%OUTPUT
%  found - True/False - found or not the command
%  curl_cmd - Path to the first found curl executable
%
    curl_cmd = '';
    
    if ispc
        [res, curl_cmd] = system('where curl');

        found = (res == 0);
        if not(found); return; end;
        
        ind = find(curl_cmd == 10, 1);
        if ~isempty(ind); curl_cmd = curl_cmd(1:ind - 1); end;

    elseif isunix
        [res, curl_path] = system('which curl');
        
        found = (res == 0);
        if not(found); return; end;

        ind = find(curl_path == 10, 1);
        if ~isempty(ind); curl_cmd = curl_path(1:ind - 1); end;
    else
        error('Unsupported OS')
    end
    return
end

function success = download_with_curl(fullurl, outdir, outfile)
% Download using the curl_command
    [found, curl_cmd] = find_curl();
    
    if not(found)
        error('Cannot download using "curl"');
    end
    
    curdir = pwd();
    cd(outdir)

    cmd = [curl_cmd, ' -O ', fullurl, ' -o ', outfile];
    res = system(cmd);
    cd (curdir);
    success = (res == 0);
end
