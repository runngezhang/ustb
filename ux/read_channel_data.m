function data = read_channel_data(filename,Nch)
% data = read_channel_data(filename,Nch)
%
% Reads channel data saved by the "channels" Texo program.
%
% Input:
% filename  - 
% Nch       - Number of channels, 32 or 64
%
% Output:
% data
if nargin < 2
    Nch = 32;
end

fp = fopen(filename,'r');
if fp == -1
    warning('could not read file! data returned is empty');
    data = [];
    return
end;
Nsamples = fread(fp,1,'int32');
    
if ~exist('data','var')
    data = zeros(Nsamples,Nch);
end    
for kk=1:Nch
    tmp = fread(fp,Nsamples,'int32');
    data(:,kk) = tmp;
    
    if kk == 32 && feof(fp)
        if Nch == 64
            warning('There is not more than 32 channels of available data');
        end            
        break;
    end
end

fclose(fp);