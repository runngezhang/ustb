function rf = correct_rf_jitter(rf,probe_id,line_density)
% rf = correct_rf_jitter(rf,probe_id,line_density)
%
% Corrects the line to line jitter observed when the line density is twice
% that of the number of elements of the probe. If the line density is equal
% to the number of elements, no correction is done.
%
% Input:
% rf                - RF signal
% probe_id          - ID of the probe
% line_density      - Line density used by the scanner (hdr.ld)
%
% Output:
% rf                - Corrected RF signal.
correct_jitter = false;
switch probe_id
    case 2    
        correct_jitter = line_density == 256;                    
    case 18    
        correct_jitter = line_density == 256;                    
    case 29
        correct_jitter = line_density == 384;
end

% Correct for 2 sample jitter
if correct_jitter
    for kk=1:2:size(rf,2)
         rf(3:end,kk,:,:) = rf(1:end-2,kk,:,:);
    end
end