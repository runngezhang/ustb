function out_dataset = stack(in_dataset)
%STACK Stacks scanlines together to produce an image

assert(isa(in_dataset(1).scan,'huff.linear_scan'),'Stack only works with LINEAR_SCAN');

out_dataset=in_dataset(1);
wb=waitbar(0,'Stacking');
x_axis=in_dataset(1).scan.x(1);
z_axis=in_dataset(1).scan.z;
for n=2:length(in_dataset)
    x_axis(n)=in_dataset(n).scan.x(1);
    out_dataset.data(:,n)=in_dataset(n).data;
end
out_dataset.scan=huff.linear_scan(x_axis.',z_axis);
close(wb);

end

