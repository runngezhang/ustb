function build_bf(use_v4sf,use_trig_tables,use_sse4,par_for)
% Function for building the planewave C++ functions for MATLAB. Several
% options are available for increase in speed.
%
% Input:
% use_v4sf          - If SSE3 extensions are available on your platform, set
% this flag to 1. Gives increase in speed but same accuracy as standard c++
% implementation.
% use_trig_tabeles  - Use lookups tables for the trigonometric functions.
% Useful for beamforming of IQ signals where the signal must phase adjusted
% accoring to the applied delay and demodulation frequency. Can give
% decreased accuracy. Currently the accuracy is set to 1.25e-4 radians.
% Increase the trig_table_sz in this file to increase accuracy.
% use_sse4          - If SSE4 is available set this flag to 1. Shouldnt
% make any difference....
% par_for           - Enable parallelization when beamforming more than one
% frame
% 
if nargin < 1
    use_v4sf = 0;
end

if nargin < 2
    use_trig_tables = 0;
end

if nargin < 3
    use_sse4 = 0;
end

if nargin < 4
    par_for = 0;
end

trig_table_sz = 50000;
currdir = pwd;
cd(fileparts(which('build_bf.m')));

try
    cfiles = ['mex' filesep 'planewave_beamforming2.cpp '];
    d_flags = '';
    ld_flags = '';
    
    sys_str = computer;
    
    if ~isempty(strfind(sys_str,'PCWIN'))
        obj_ext = 'obj';
        d_flags = [d_flags '-D_WIN_ '];
        if par_for
            disp('# Parallelization for beamforming of multiple frames is enabled.')
            d_flags = [d_flags '-D_PAR_FOR_ '];
        end
    elseif ~isempty(strfind(sys_str,'GLNX'))        
        if par_for
            disp('# Parallelization for beamforming of multiple frames is enabled.')
            d_flags = [d_flags '-D_PAR_FOR_ '];
            ld_flags = 'LDOPTIMFLAGS=''-O -fopenmp'' ';
        end
        
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
        
        if par_for
            optim_flags = [optim_flags '  -fopenmp'];
        end
        
        d_flags = [d_flags sprintf('-D_UNIX_ CXXOPTIMFLAGS=''%s'' ',optim_flags)];
        
        obj_ext = 'o';
    else
        error('Not supported on MAC systems');
    end
    
    if use_sse4 && use_v4sf
        disp('# Make sure you CPU supports SSE4');
        d_flags = [d_flags '-D_HAS_SSE4_ '];
    end      
    
    if use_trig_tables
        disp('# Using trig tables for sin and cos')
        d_flags = [d_flags '-D_USE_TRIG_TABLES_ -' sprintf('D_TRIG_TABLE_SZ_=%d ',trig_table_sz)];
    end
    
    if use_v4sf
        disp('# Make sure you CPU supports SSE3');
        d_flags = [d_flags '-D_USE_V4SF_ '];
        
        eval(['mex -c ' d_flags ' mex' filesep 'planewave_bf_v4sf.cpp']);
        eval(['mex -c ' d_flags ' mex' filesep 'planewave_beamforming2.cpp']);
        
        cfiles = sprintf('planewave_beamforming2.%s planewave_bf_v4sf.%s',obj_ext,obj_ext);
    else                
        eval(['mex -c ' d_flags ' mex' filesep 'planewave_bf.cpp']);
        eval(['mex -c ' d_flags ' mex' filesep 'planewave_beamforming2.cpp']);
        
        cfiles = sprintf('planewave_beamforming2.%s planewave_bf.%s',obj_ext,obj_ext);
    end
       
    mex_str = ['mex ' ld_flags cfiles];
    eval(mex_str);
catch err
    cd(currdir);
    rethrow(err);
end
delete('*.o');
cd(currdir);
    

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