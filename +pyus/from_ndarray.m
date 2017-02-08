function res = from_ndarray(data)
    keySet = {'uint8','float32','float64','complex64','complex128'};
    valueSet = {'uint8','single','double','complex','complex'};
    
    mapper = containers.Map(keySet,valueSet);
    type = mapper(char(data.dtype.name));

    shapedata = typecast(uint8(py.memoryview(py.numpy.array(data.shape)).tobytes),'uint64');
    
    if ~strcmp(type,'complex')
        bytedata = uint8(py.memoryview(py.numpy.ascontiguousarray(data).T).tobytes);
        bytedata = typecast(bytedata,type);
        res = reshape(bytedata,shapedata);
    else
        realdata = uint8(py.memoryview(py.numpy.ascontiguousarray(data.real.T,'float64')).tobytes);
        realdata = typecast(realdata,'double');
        imagdata = uint8(py.memoryview(py.numpy.ascontiguousarray(data.imag.T,'float64')).tobytes);
        imagdata = typecast(imagdata,'double');
        bytedata = complex(realdata,imagdata);
        res = reshape(bytedata,shapedata);
    end

end

