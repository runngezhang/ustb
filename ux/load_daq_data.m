function [hdr, RFframe] = load_daq_data(data_path, chanls, frames, reRoute)
% Read the specified channel data from .daq files
% Usage:
% [hdr, RFframe] = load_daq_data(data_path, chanls, frames,reRoute);
%       data_path:  strings showing the path to the data folder
%       chanls:     chanels of interest
%       frames:     index of frames of interest, 1,2...
%       reRoute:    true: transducer element, false: DAQ element
% Example:
%       Reading single frame
%       [header, RFframe] = readRFframe('c:\foldername\', ones(1,128), 2, true);
%       will return the second frame as a 2D matrix in RFframe
% 
% Author: Reza Zahiri Azar, March 2008, University of British Columbia
% Changed by Thor Andreas Tangen, January 2010, NTNU


% Signal assignement to transducer element
if (reRoute)
    chRout = [0 16 32 48 64 80 96 112 1 17 33 49 65 81 97 113 2 18 34 50 66 82 98 114 ...
              3 19 35 51 67 83 99 115 4 20 36 52 68 84 100 116 5 21 37 53 69 85 101 117 ...
              6 22 38 54 70 86 102 118 7 23 39 55 71 87 103 119 8 24 40 56 72 88 104 120 ...
              9 25 41 57 73 89 105 121 10 26 42 58 74 90 106 122 11 27 43 59 75 91 107 ...
              123 12 28 44 60 76 92 108 124 13 29 45 61 77 93 109 125 14 30 46 62 78 94 ...
              110 126 15 31 47 63 79 95 111 127];
else    % show channel signals
    chRout = 0:127;  
end

% open channel file
            
for ii = 1:128        
    if (chanls(1 + chRout(ii)) == 1)
        % create the file name
        chlInd = ii-1;
        if (chlInd<10)
            tag = ['00',num2str(chlInd)];
        elseif (chlInd<100)
            tag = ['0',num2str(chlInd)];
        else
            tag = num2str(chlInd);
        end
        
        filename = fullfile(data_path,['CH',tag,'-0000.daq']);
        fid = fopen(filename,'r');
        
        if fid < 0            
            error('Could not load file %s',filename);
        end
        
        % read header
        hdr = fread(fid, 2, 'int32');
        numFrame = hdr(1);      % number of frames acquired
        lLength  = hdr(2);      % length of each line in samples
        
        for kk=1:length(frames)                        
            if frames(kk) > numFrame
                warning('Specified frame number exceeds the number of acquired frames');
                break;
            end
            
            if ~exist('RFframe','var')
                RFframe = zeros(lLength,128,length(frames));
            end
            
            if kk == 1
                frame_skip = frames(kk);
            else
                frame_skip = frames(kk) - frames(kk-1);
            end                        

            if fseek( fid, (frame_skip - 1) * ( lLength * 2), 'cof') < 0                               
                error(ferror(fid))
            end
            
            % read channel data and correct the mapping
            dataformat = 'int16';
            ind = 1 + chRout(ii);
            RFframe(:,ind,kk) = fread(fid, lLength, dataformat );   % the actual data                       
            fclose(fid);
        end;
    end;
    % close channel file    
end

