function val = getAvailableMemory
if isunix
    [~, val] = system('grep MemAvailable /proc/meminfo');
    val = sscanf(val, 'MemAvailable:   %d kB');
elseif ispc
    [~, sys] = memory;
    val = sys.PhysicalMemory.Available;
else
    disp('Platform not supported')
end
end