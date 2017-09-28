function check_memory(bytes)
%CHECK_MEMORY Checks if there is room for more data in the RAM

[user,sys] = memory;

safety_factor = 0.3;

assert(sys.PhysicalMemory.Available>(1 + safety_factor)*bytes,'Not available RAM for the new data. Aborting process.');

end
