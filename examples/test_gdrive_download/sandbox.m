clear all
close all
clc

id = '19OyvPCP4qUiTECFpUe8r3r_Ys2281j0N'; % unique file id on google drive
url = 'https://drive.google.com';

% We have to send a first GET request to download the "warning download page"
response = send(matlab.net.http.RequestMessage(), strcat(url, '/uc', '?', 'export=download', '&id=', id));

% Then we extract the cookies from the response
setCookie = response.getFields('Set-Cookie');
cookieInfo = setCookie.convert();

% We prepare the second GET request, which starts the file download
request = matlab.net.http.RequestMessage();

% We scan the received cookies
for cookie = [cookieInfo.Cookie]
    
    % We look for the cookie whose fields starts with "download warning"
    if startsWith(cookie.Name, 'download_warning')
        key = cookie.Value;
        request = addFields(request, 'Cookie', matlab.net.http.Cookie(cookie.Name, cookie.Value));
    end
end

% We prepare a http FileConsumer object to take care of the download
consumer = matlab.net.http.io.FileConsumer('Verasonics_P2-4_apical_four_chamber_subject_1.uff');

% We specify a http ProgressMonitor object to display the download progress 
% on the command window

opts = matlab.net.http.HTTPOptions('ProgressMonitorFcn', @progressMonitor, 'UseProgressMonitor',true);
send(request, strcat(url, '/uc', '?', 'export=download', '&id=', id, '&confirm=', key), opts, consumer);