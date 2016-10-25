function res = fromNdarray(data)
    keySet = {'uint8','float32','float64','complex64','complex128'};
    valueSet = {'uint8','single','double','complex','complex'};
    mapper = containers.Map(keySet,valueSet);

    type = mapper(char(data.dtype.name));

    shapedata = typecast(uint8(py.memoryview(py.numpy.array(data.shape)).tobytes),'uint64');
    tmp = shapedata(1);
    shapedata(1) = shapedata(2);
    shapedata(2) = tmp;
    
    if ~strcmp(type,'complex')
        bytedata = uint8(py.memoryview(py.numpy.ascontiguousarray(data)).tobytes);
        bytedata = typecast(bytedata,type);
        res = reshape(bytedata,shapedata);
    else
        realdata = uint8(py.memoryview(py.numpy.ascontiguousarray(data.real,'float64')).tobytes);
        realdata = typecast(realdata,'double');
        imagdata = uint8(py.memoryview(py.numpy.ascontiguousarray(data.imag,'float64')).tobytes);
        imagdata = typecast(imagdata,'double');
        bytedata = complex(realdata,imagdata);
        res = reshape(bytedata,shapedata);
    end

end

