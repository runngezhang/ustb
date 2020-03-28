clear all
close all
clc

file_id = 2;    % 1. large file. requires download confirmation
                % 2. small file. does not require download confirmation

switch file_id
    case 1
        id = '19OyvPCP4qUiTECFpUe8r3r_Ys2281j0N'; % file id on google drive
        name = 'Verasonics_P2-4_apical_four_chamber_subject_1.uff';
    case 2
        id = '1LDu9WGeYvOFY--TnIHrUHBsvranQsPKf';
        name = 'PICMUS_carotid_long.uff';
end

url = 'https://drive.google.com';
opts = matlab.net.http.HTTPOptions('ProgressMonitorFcn', @progressMonitor, ...
    'UseProgressMonitor',true);

% We send a first GET request
response = send(matlab.net.http.RequestMessage(), strcat(url, '/uc', '?', ...
    'export=download', '&id=', id), opts);

if response.StatusCode == 200
    % First, we check that the response from the server was OK
    
    if strcmp(response.Body.ContentType.Type, 'application')
        % If the content of the first response is of type 'application-
        % octet-stream', it means that the file was small enough not to 
        % prompt the warning download message. Therefore, we have already
        % downloaded the file. We just need to save it.
        
        fid = fopen(fullfile(ustb_path, name), 'w');
        fwrite(fid, response.Body.Data);
        fclose(fid);
        
    elseif strcmp(response.Body.ContentType.Type, 'text')
        % If the content of the first response is of type 'text-html', it
        % means that the file was large enough to prompt the warning
        % download message. Therefore, we need to send a confirm request to
        % start the download
        
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
        
        % We prepare a http FileConsumer object to take care of writing the 
        % file on disk
        consumer = matlab.net.http.io.FileConsumer(fullfile(ustb_path, name));
        
        % We send the second GET request and start the file download
        response = send(request, strcat(url, '/uc', '?', 'export=download', ...
            '&id=', id, '&confirm=', key), opts, consumer);
    else
        error(1, 'The response from the server was not recognised');
    end

else
    error(response.StatusCode, 'The HTTP request failed with error %d\n', ...
        response.StatusCode);
end