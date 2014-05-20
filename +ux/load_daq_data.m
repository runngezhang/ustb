function varargout = load_daq_data(datapath,varargin)

% TODO: Allow for number of elements other than 128!
Nch = 128;

p = inputParser();
p.CaseSensitive = false;
p.addParameter('Reroute',false);
% p.addParameter('Size',[]);
p.addParameter('WriteToFile',[]);
p.addParameter('Precision','int16');
p.addParameter('Channels',0:Nch-1);
p.addParameter('FileNo',0);
p.addParameter('FileFormat','6.0');
p.addParameter('Frames',[]);
p.parse(varargin{:});
res = p.Results;

for fileno = res.FileNo(:)',
    ch = res.Channels;
    if str2double(res.FileFormat) < 6.0,
        daqfilename0 = fullfile(datapath,sprintf('CH%03d-%04d.daq',ch(1),fileno));
        header_length = 2; % bytes
        hdr = read_daq_header(daqfilename0,header_length);
    else
        daqfilename0 = fullfile(datapath,sprintf('CH%03d.daq',ch(1)));
        header_length = 3; % bytes
        hdr = read_daq_header(daqfilename0,header_length);
    end
    if isempty(res.Frames),
        res.Frames = 1:hdr(1);
    end
    Nrf = hdr(2);
    
    if true,%isempty(fsize),
        fsize = [Nrf,numel(ch),numel(res.Frames)];
    end
    
    % Signal assignement to transducer element
    chRout = 0:Nch-1;
    if res.Reroute,
        Q = 16;
        chRout = reshape(reshape(chRout',Q,[])',1,[]);
    end
    
    if res.WriteToFile,
        if numel(res.FileNo) > 1,
            [outdir,fname,ext] = fileparts(res.WriteToFile);
            fmtstr = sprintf('%%s%%0%dd%%s',floor(log10(max(res.FileNo)))+1);
            matfilename = sprintf(fmtstr,fullfile(outdir,fname),fileno,ext);
        else
            matfilename = res.WriteToFile;
        end
        S = matfile(matfilename,'Writable',true);
    else
        S = struct();
    end
    S.RFframe = zeros(fsize,res.Precision);
    
    [frames_,~,frind] = unique(res.Frames(:));
    frame_density = numel(frames_)/(frames_(end)-frames_(1)+1);
    
    % open channel file
    for jj = 1:numel(ch),
        % filename = fullfile(datapath,sprintf('CH%03d-%04d.daq',ch(jj),res.FileNo));
        filename = fullfile(datapath,sprintf('CH%03d.daq',ch(jj)));
        fid = fopen(filename,'r');
        hcl = onCleanup(@()fclose(fid));
        
        if frame_density > 1/20,
            % Most efficient to read chunk of frames, and sort them out later
            frame_skip = header_length*4 + 2*Nrf*(frames_(1)-1); % skip header too, 8 bytes
            fseek(fid,frame_skip,'bof');
            T = fread(fid,[Nrf,frames_(end)-frames_(1)+1],['*',res.Precision]);
            S.RFframe(1:Nrf,1+chRout(1+ch(jj)),1:numel(frind)) = permute(T(:,frind),[1,3,2]);
        else
            % Most efficient to loop
            fseek(fid,header_length*4,'bof');
            frame_skip = 2*Nrf*(diff([0;frames_])-1);
            for kk = 1:numel(frames_),
                fseek(fid,frame_skip(kk),'cof');
                S.RFframe(1:Nrf,1+chRout(1+ch(jj)),kk) = fread(fid,[Nrf,1],['*',precision]);
            end
        end
    end 
end

if nargout,
    varargout{1} = S.RFframe;
end

    function h = read_daq_header(filename,hdrlen)
        fid0 = fopen(filename,'r');
        h = fread(fid0, hdrlen, 'int32');
        h = h(end-1:end);
        fclose(fid0);
    end

end
