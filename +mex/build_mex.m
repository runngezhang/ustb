function a=build_mex(filename)

use_sse4 = true;
use_v4sf = false;

try
    d_flags = '';
    ld_flags = '';
    
    sys_str = computer;
    
    if ~isempty(strfind(sys_str,'PCWIN'))
        obj_ext = 'obj';
        d_flags = [d_flags '-D_WIN_ '];
        
    elseif ~isempty(strfind(sys_str,'GLNX'))        
        ld_flags = 'LDOPTIMFLAGS=''-O -fopenmp'' ';
        
        [main_ver, minor_ver, rev_ver] = get_gcc_version();
        
        if main_ver < 4
            error('Must have gcc version 4 or higher');
        end
        
        if main_ver == 4 && minor_ver < 3 && use_sse4
            warning('Must have gcc version 4.3 or higher in order to use SSE4. Disabling SSE4');
            use_sse4 = 0;
        end
                
        optim_flags = '-O2 -DNDEBUG -msse -msse2';
        
        if use_sse4
            optim_flags = [optim_flags ' -msse4.1'];
        end
        
        optim_flags = [optim_flags '  -fopenmp'];
        
        d_flags = [d_flags sprintf('-D_UNIX_ CXXOPTIMFLAGS=''%s'' ',optim_flags)];
        
        obj_ext = 'o';
    else
        error('Not supported on MAC systems');
    end
    
    if use_sse4 && use_v4sf
        disp('# Make sure you CPU supports SSE4');
        d_flags = [d_flags '-D_HAS_SSE4_ '];
    end      
    
    eval(['mex -c ' d_flags ' source' filesep filename '.cpp']);
        
    cfiles = sprintf('%s.%s',filename,obj_ext);
       
    mex_str = ['mex ' ld_flags cfiles];
    eval(mex_str);
catch err
    %cd(currdir);
    rethrow(err);
end
delete('*.o');
    

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
    break;
end