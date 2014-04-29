function [data,fs,type,probe_id] = load_signal(filename)
% [data,fs,type,probe_id] = load_signal(filename)
%
% Output parameters
% data          - Data, 4 dimensional array: [Length,LineDensity,EnsembleSize,NumFrames]
% fs            - Sampling frequency [Hz]
% type          - Data type. See ulterius_def.h.
% probe_id      - Id of the probe used

[x header] = load_ux_signal(filename);

line_density = header.ld;
ensemble_size = header.w/line_density;
num_frames = size(x,3);

data = zeros(header.h,line_density,ensemble_size,num_frames);

for k=1:num_frames
    for n=1:ensemble_size
        data(:,:,n,k) = x(:,n:ensemble_size:end,k);
    end
end

fs = header.sf;
type = header.type;
probe_id = header.probe;

function [x header frame_numbers data_type] = load_ux_signal(filename)

[path,name,ext] = fileparts(filename);

fp = fopen(filename,'r');
hdr_fields = ux_header();

if fp > 0 
    header_values = double(fread(fp,length(hdr_fields),'int32'));
    if length(header_values) < length(hdr_fields)
        error('Error reading file');        
    end
    
    for k=1:(length(hdr_fields))
        header.(hdr_fields{k}) = header_values(k);
    end    
    
    precision = get_precision(header.type);
    data_type = get_datatype(header.type);
    
    num_samples = header.w*header.h;
    has_appended_data = 0;
    has_prepended_data = 0;
    
    if strcmp(ext,'.bpr') || strcmp(ext,'.rf')
        has_appended_data = 1;
    end
    header.frames = max(header.frames,1);
    
    x = zeros(header.h,header.w,header.frames);
    frame_numbers = -1*ones(1,header.frames);
    all_samples = fread(fp,inf,precision);
    
    if length(all_samples) ~= header.frames*num_samples
        has_appended_data = 1;
    else
        has_appended_data = 0;
    end

    frame_start = 0;
    for k=1:header.frames
        
        if has_prepended_data
            frame_numbers(k) = fread(fp,1,'int32');
        end
        
        tmp = all_samples(frame_start + (1:num_samples));                        
        frame_start = frame_start + num_samples;
        
        if has_appended_data
            if strcmp(precision,'int16')
                frame_start = frame_start + 2;
            elseif strcmp(precision,'uchar')
                frame_start = frame_start + 4;
            elseif strcmp(precision,'float32')
                frame_start = frame_start + 1;
            end
        end
        
        tmp = reshape(tmp,header.h,header.w);    
                
        x(:,:,k) = tmp;
    end
    fclose(fp);
end

function hdr = ux_header()

% data type - data types can also be determined by file extensions
hdr{1} = 'type';
% number of frames in file
hdr{2} = 'frames';	
% width - number of vectors for raw data, image width for processed data    
hdr{3} = 'w';
% height - number of samples for raw data, image height for processed data
hdr{4} = 'h';
% data sample size in bits    
hdr{5} = 'ss';
% roi - upper left (x)
hdr{6} = 'ulx';
% roi - upper left (y)
hdr{7} = 'uly';
% roi - upper right (x)
hdr{8} = 'urx';
% roi - upper right (y)
hdr{9} = 'ury';
% roi - bottom right (x)
hdr{10} = 'brx';
% roi - bottom right (y)
hdr{11} = 'bry';
% roi - bottom left (x)
hdr{12} = 'blx';
% roi - bottom left (y)
hdr{13} = 'bly';
% probe identifier - additional probe information can be found using this id
hdr{14} = 'probe';
% transmit frequency
hdr{15} = 'txf';
% sampling frequency
hdr{16} = 'sf';    
% data rate - frame rate or pulse repetition period in Doppler modes
hdr{17} = 'dr';
% line density - can be used to calculate element spacing if pitch and native # elements is known
hdr{18} = 'ld';
% extra information - ensemble for color RF
hdr{19} = 'extra'; 

function precision = get_precision(type)

switch(type)
    case 1
        precision = 'uchar';
    case 2
        precision = 'uchar';
    case 4
        precision = 'uchar';
    case 8
        precision = 'float32';
    case 16
        precision = 'int16'; 
    otherwise
        precision = 'uchar';
end

function data_type = get_datatype(type)

switch(type)
    case 1
        data_type = 'screen';
    case 2
        data_type = 'b-pre';
    case 4
        data_type = 'b-post';
    case 8
        data_type = 'b-post32';
    case 16
        data_type = 'rf';
    otherwise
        data_type = '';
end