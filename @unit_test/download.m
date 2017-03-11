function download(h, filename, url, local_path )
%DOWNLOAD Downloads a file from an URL to a local path (if needed)

    % check if the path is there 
    a=dir([local_path]);
    if not(numel(a))
        % if not create it
        mkdir(local_path);
    end
    
    % check if file is in local path
    a=dir([local_path filename]);
    if not(numel(a))
        % if not download it
        disp(['Downloading ' url filename '. This may take a while.']);
        websave([local_path filename],[url filename]);
    end

end

