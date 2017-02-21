%% Using HUFF class to read
%
% This example shows how to interact with the huff class to read data from a hdf5 
% file according to the specification HUFF v0.0.1.

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2015/02/03 $

%% create a huff object
recording=huff('recording.h5','r');

%% load everything inside the file
recording.read();

%% plot all spatial reconstructions
for n=1:length(recording.spatial_reconstruction)
    recording.spatial_reconstruction{n}.show();
end

%% reconstruct all images from the datasets
recons=reconstruction('New reconstruction',recording.spatial_reconstruction{1});    % new reconstruction copied from previous one
for n=1:length(recording.spatial_reconstruction)
    recording.ultrasound_signal{n}.image_reconstruction(recons);
    recons.show();
end
