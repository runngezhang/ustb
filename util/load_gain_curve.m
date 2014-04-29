function [gain_curve,z_mm] = load_gain_curve(probe_id,application_id,c,fs,N,fname)
% gain_curve = load_gain_curve(probe_id,application_id,dz_mm,N,fname)
%
% Loads the gain curve for teh specified probe and application
%
% Input:
% probe_id          - ID of probe
% application_id    - ID of application, e.g. SURF_TCI, B_MODE
% c                 - Speed of sound
% fs                - Sampling frequency
% N                 - Number of samples
% fname             - XML where the gain curves are given (Optional)
% 
% Output:
% gain_curve        - The gain curve
% z_mm              - Depth axis in mm

if nargin < 6    
    ustb_path = fileparts(which('ustb'));
    fname = fullfile(ustb_path,'config_files','gain_curves.xml');
end

doc = xmlread(fname);
node = doc.getDocumentElement;
all_probes = node.getElementsByTagName('probe');

pp = [];
for k = 0:all_probes.getLength-1
    probe_node = all_probes.item(k);
    id = str2double(char(probe_node.getAttribute('id')));
    
    if probe_id == id
        app_nodes = probe_node.getElementsByTagName('application');

        for n = 0:app_nodes.getLength-1
            app_node = app_nodes.item(n);
            id = str2double(char(app_node.getAttribute('id')));
            
            if application_id == id
                
                gain_nodes = app_node.getElementsByTagName('gaincurve');
                for l = 0:gain_nodes.getLength-1
                    gain_node = gain_nodes.item(l);
                    id = str2double(char(gain_node.getAttribute('id')));
                    if id == 0
                        
                        pp = str2num(char(gain_node.getAttribute('coeffs')));
                        zlims = str2num(char(gain_node.getAttribute('zlims')));
                        
                    end
                end
            end
        end
    end
end

dz_mm = 1e3*c/fs/2;
z_mm = (0:(N-1))'*dz_mm;
z_ix = find_indx(z_mm,zlims);
    
if isempty(pp)
    gain_curve = ones(size(z_mm));
else    
    gain_curve = polyval(pp,z_mm);
    gain_curve = 1./gain_curve;

    dg = gain_curve(z_ix(1))./z_ix(1);
    gain_curve(1:z_ix(1)) = dg.*(0:(z_ix(1)-1));
end