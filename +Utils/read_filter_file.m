function [filters,header] = read_filter_file(filename)


fp = fopen(filename,'r');

for k=1:6
    line = fgetl(fp);
    header{k} = line;    
end

filter_cnt = 1;
line_cnt = 1;
while 1
    line = fgetl(fp);

    if ~ischar(line)
        break;
    end
    
    switch line_cnt
        case 1
            filter.id = sscanf(line,'<%d>');
        case 2
            filter.type = line;
        case 3
            filter.desc = line;
        case 4
            filter.num_taps = sscanf(line,'%d');
        case 5
            filter.taps = sscanf(line,'%d');
        case 6
            filter.scale = sscanf(line,'%d');
    end
    
    line_cnt = line_cnt + 1;
    
    if line_cnt == 7
        filter_cnt = filter_cnt + 1;
        line_cnt = 1;
        filters{filter.id+1} = filter;
    end
end

fclose(fp);