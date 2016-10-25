function setup_pyus(path_to_pyus)

if count(py.sys.path, path_to_pyus) == 0
    insert(py.sys.path,int32(0),path_to_pyus);
end

end

