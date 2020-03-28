function download(file, url)
%DOWNLOAD Downloads a file from an URL to a local path (if needed)
%
% NOTE: compatible with

[path, name, ext] = fileparts(file);

% Check that the file has not been downloaded previously
if not(exist(file,  'file'))
    fprintf(1, 'Downloading %s%s\nfrom %s\nto %s\n', name, ext, url, path)
    
    % Create folder if it does not exist
    if ~exist(path, 'dir')
        mkdir(path)
    end
    
    % Prepare a HTTP option object were we specify to use a custom progress
    % monitor, which updates the user on the amount of downloaded data
    opts = matlab.net.http.HTTPOptions('ProgressMonitorFcn', ...
        @tools.progressMonitor, 'UseProgressMonitor',true);
    
    % We send a first GET request
    response = send(matlab.net.http.RequestMessage(), url, opts);
    
    % First, we check that the response from the server was OK
    if response.StatusCode == 200
        
        % If the content of the first response is of type 'application-
        % octet-stream' or is not specified, it means that we have already
        % downloaded the file in the first request. NOTE: the order in
        % which the two clauses are checked is important
        if isempty(response.Body.ContentType) || ...
                strcmp(response.Body.ContentType.Type, 'application')
            
            % We just need to save the file
            fid = fopen(file, 'w');
            fwrite(fid, response.Body.Data);
            fclose(fid);
            
        % If the content of the first response is of type 'text-html', it
        % means that the file was large enough to prompt the warning
        % download message. Therefore, we need to send a confirm request to
        % start the download
        elseif strcmp(response.Body.ContentType.Type, 'text')        
            
            % First, we prepare a second GET request, which will start the file
            % download
            request = matlab.net.http.RequestMessage();
            
            % Then we extract the cookies from the response
            setCookie = response.getFields('Set-Cookie');
            cookieInfo = setCookie.convert();
            
            % We look for the cookie whose field starts with "download warning"
            for cookie = [cookieInfo.Cookie]
                
                if startsWith(cookie.Name, 'download_warning')
                    
                    key = cookie.Value;
                    request = addFields(request, 'Cookie', ...
                        matlab.net.http.Cookie(cookie.Name, cookie.Value));
                end
            end
            
            % We send the second GET request and start the file download
            response = send(request, strcat(url, '&confirm=', key), opts);
            
            % We save the file
            fid = fopen(file, 'w');
            fwrite(fid, response.Body.Data);
            fclose(fid);
        else
            error('Unknown content type!');
        end
        
    else
        error('The HTTP request failed with error %d', response.StatusCode);
    end
end
end