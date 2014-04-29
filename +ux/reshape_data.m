function rf_out = reshape_data(rf_in,zone_lengths)
% rf_out = reshape_data(rf_in,focus_samples)
%
% Function for restacking data in the case of multiple foci and ensemble
% size greater than 1.
%
% Input:
% rf_in         - RF data to be restacked. Has dimensions [Length LineDensity EnsembleSize NumFrames]
% zone_lengths - Length of each focal zone in number of samples . Vector
% of length NumFoci
%
% rf_out        - Restacked data

rf_out = zeros(size(rf_in));
dims = size(rf_in);

ensemble_size = dims(3);
num_frames = dims(4);
num_foci = length(zone_lengths);

focus_starts = [0 cumsum(zone_lengths(1:end-1))];

indx1 = [];
indx2 = [];    
for e=1:ensemble_size        
    for f=1:num_foci            
        start_indx1 = focus_starts(f) + (e-1)*dims(1);                        
        start_indx2 = ensemble_size*focus_starts(f) + (e-1)*zone_lengths(f);                  
        indx1 = [indx1 start_indx1 + (1:zone_lengths(f))];            
        indx2 = [indx2 start_indx2 + (1:zone_lengths(f))];            
    end
end

offset = ones(length(indx1),1)*2*(0:(dims(2)-1))*dims(1);
indx1 = offset(:)' + repmat(indx1,1,dims(2));
indx2 = offset(:)' + repmat(indx2,1,dims(2));

rf1 = zeros(dims(1),dims(2)*dims(3));
for n=1:num_frames    
        
    for m=1:ensemble_size
        rf1(:,m:ensemble_size:end) = rf_in(:,:,m,n);
    end                   
       
    rf1(indx1) = rf1(indx2);           
    
    for m=1:ensemble_size
        rf_out(:,:,m,n) = rf1(:,m:ensemble_size:end);
    end    
end