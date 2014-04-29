function [x header params actual_frames] = load_ux_signal(filename,num_frames,start_frame)

[path,name,ext] = fileparts(filename);
xml_file = [path filesep name '.xml'];
if exist(xml_file,'file') && nargout == 3
    params = params_from_xml(xml_file);
else
    params = [];
end

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
    
    actual_frames = header.frames;
    
    if header.type == 16
        ensemble = ceil(header.w/header.ld);
        header.w = header.w/ensemble;
    elseif header.type == 512
        ensemble = header.extra;
    else
        ensemble = 1;
    end        
    
    [precision,bytes_per_sample,col_major] = get_precision(header.type);
    data_type = get_datatype(header.type);    
    
    num_samples = header.w*header.h*ensemble;
    has_appended_data = 0;
    has_prepended_data = 0;
    
    if strcmp(ext,'.bpr') || strcmp(ext,'.rf')
        has_appended_data = 1;
    end
    
    header.frames = max(header.frames,1);    
    if nargin < 3
        start_frame = 1;        
    end
    
    header.frames = header.frames - (start_frame-1);
    
    if nargin >= 2
        header.frames = min(header.frames,num_frames);    
    end
    
    x = zeros(header.h,header.w,ensemble,header.frames,precision);
    frame_numbers = ones(1,header.frames);
        
    % Seek to the reading position
    seek_offset = (start_frame-1)*(num_samples*bytes_per_sample + has_appended_data*4 + has_prepended_data*4);
    fseek(fp,seek_offset,'cof');        
    
    for k=1:header.frames          
        if has_prepended_data
            frame_numbers(k) = fread(fp,1,'int32') + start_frame - 1;
        end
               
        tmp = fread(fp,num_samples,precision);

        if length(tmp) < num_samples
            error('Error reading file');
        end
        
        if has_appended_data
            %frame_numbers(k) = fread(fp,1,'int32');
        end
        
        if col_major
            tmp = reshape(tmp,header.h,header.w*ensemble);            
        else            
            tmp = reshape(tmp,header.w,header.h*ensemble);                                    
            tmp = tmp';
        end
        
        for n=1:ensemble
            x(:,:,n,k) = tmp(:,n:ensemble:end);
        end
    end
    fclose(fp);
else
    error('ERROR: Could not open file %s',filename);
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

function [precision,bytes,col_major] = get_precision(type)

col_major = true;
switch(type)
    case 1
        precision = 'uchar';
        bytes = 1;
    case 2
        precision = 'uchar';
        bytes = 1;
    case 4
        precision = 'uchar';
        col_major = false;
        bytes = 1;
    case 8
        precision = 'int32';
        bytes = 4;
    case 16
        precision = 'int16'; % 16 bit RF data
        bytes = 2;
    case 128
        precision = 'uint16'; % 16 bit envelope
        bytes = 2;
    case 256
        precision = 'float32';
        bytes = 4;
    case 512 % Color RF
        precision = 'int16';
        bytes = 2;    
    case 16384
        precision = 'uchar';
        col_major = false;
        bytes = 1;
    otherwise
        precision = 'uchar';
        bytes = 1;
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
    case 128
        data_type = 'env';
    case 256
        data_type = 'env';
    case 512
        data_type = 'crf';    
    case 16384
        data_type = 'elasto_over';
    otherwise
        data_type = '';
end