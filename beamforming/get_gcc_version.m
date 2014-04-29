function [main_ver, minor_ver, rev_ver] = get_gcc_version()

main_ver = 0;
minor_ver = 0; 
rev_ver = 0;
[ret str] = system('gcc --version');
if ret == 127
    error('gcc not present on system');
end

remain = str;
while true 
    [tok remain] = strtok(remain,' ');
    if isempty(tok)
        break;
    end
    ver_nums = sscanf(tok,'%d.%d.%d');
    if isempty(ver_nums) || length(ver_nums) ~= 3
        continue
    end
    
    main_ver = ver_nums(1);
    minor_ver = ver_nums(2);
    rev_ver = ver_nums(3);
end