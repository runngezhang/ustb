function probe_info = get_probe_info(probe_id,probes_lst_path,list_probes)
% Function for retrieving probe information from the probes.lst file.
%
% Input:
% probe_id          - ID of probe
% probes_lst_path   - Where to find the probes.lst file

if nargin < 2
    probes_lst_path = ''; % look on matlab path
end

filename = fullfile(probes_lst_path,'probes.lst');

if nargin  < 3
    list_probes = false;
end

probe_info = [];

fp = fopen(filename,'r');

if fp==-1
    error('Could not find probe information file');    
end

probe_fields = get_probe_fields();

while ~feof(fp)    
    for k=1:length(probe_fields)        
        probe_info.(probe_fields{k}.name) = fread(fp,probe_fields{k}.count,probe_fields{k}.precision);
    end    
    if list_probes
        fprintf('Name : %s, ID : %d, Type : %d\n',probe_info.name,probe_info.id,probe_info.type);
    end
    if probe_id == probe_info.id        
        break;
    end        
end

int8_name = int8(probe_info.name);
indx = find(int8_name==0,1,'first');
probe_info.name = probe_info.name(1:(indx-1))';

fclose(fp);

function probe_fields = get_probe_fields()

probe_fields{1}.name = 'name';
probe_fields{1}.precision = 'uint8=>char';
probe_fields{1}.count = 80;

probe_fields{2}.name = 'id';
probe_fields{2}.precision = 'int32';
probe_fields{2}.count = 1;

probe_fields{3}.name = 'type';
probe_fields{3}.precision = 'int32';
probe_fields{3}.count = 1;

probe_fields{4}.name = 'elements';
probe_fields{4}.precision = 'int32';
probe_fields{4}.count = 1;

probe_fields{5}.name = 'pitch';
probe_fields{5}.precision = 'int32';
probe_fields{5}.count = 1;

probe_fields{6}.name = 'radius';
probe_fields{6}.precision = 'int32';
probe_fields{6}.count = 1;

probe_fields{7}.name = 'freq';
probe_fields{7}.precision = 'int32';
probe_fields{7}.count = 1;

probe_fields{8}.name = 'freqBW';
probe_fields{8}.precision = 'int32';
probe_fields{8}.count = 1;

probe_fields{9}.name = 'transmitoffset';
probe_fields{9}.precision = 'int32';
probe_fields{9}.count = 1;

probe_fields{10}.name = 'maxsteerangle';
probe_fields{10}.precision = 'int32';
probe_fields{10}.count = 1;

probe_fields{11}.name = 'maxfocusdistance';
probe_fields{11}.precision = 'int32';
probe_fields{11}.count = 1;

probe_fields{12}.name = 'minlineduration';
probe_fields{12}.precision = 'int32';
probe_fields{12}.count = 1;

probe_fields{13}.name = 'min_focus_distance_doppler';
probe_fields{13}.precision = 'int32';
probe_fields{13}.count = 1;

probe_fields{14}.name = 'pin_offset';
probe_fields{14}.precision = 'uchar';
probe_fields{14}.count = 1;

probe_fields{15}.name = 'unusedC';
probe_fields{15}.precision = 'schar';
probe_fields{15}.count = 1;

probe_fields{16}.name = 'unusedS';
probe_fields{16}.precision = 'int16';
probe_fields{16}.count = 1;

probe_fields{17}.name = 'motorFOV';
probe_fields{17}.precision = 'int32';
probe_fields{17}.count = 1;

probe_fields{18}.name = 'biospy_depth';
probe_fields{18}.precision = 'int32';
probe_fields{18}.count = 1;

probe_fields{19}.name = 'biopsy_angle';
probe_fields{19}.precision = 'int32';
probe_fields{19}.count = 1;

probe_fields{20}.name = 'biopsy_distance';
probe_fields{20}.precision = 'int32';
probe_fields{20}.count = 1;

probe_fields{21}.name = 'options';
probe_fields{21}.precision = 'int32';
probe_fields{21}.count = 1;

probe_fields{22}.name = 'motor_steps';
probe_fields{22}.precision = 'int32';
probe_fields{22}.count = 1;

probe_fields{23}.name = 'motor_radius';
probe_fields{23}.precision = 'int32';
probe_fields{23}.count = 1;

probe_fields{24}.name = 'motor_min_time_between_pulses';
probe_fields{24}.precision = 'uint16';
probe_fields{24}.count = 1;

probe_fields{25}.name = 'motor_home_method';
probe_fields{25}.precision = 'schar';
probe_fields{25}.count = 1;

probe_fields{26}.name = 'biopsy_width';
probe_fields{26}.precision = 'uchar';
probe_fields{26}.count = 1;
