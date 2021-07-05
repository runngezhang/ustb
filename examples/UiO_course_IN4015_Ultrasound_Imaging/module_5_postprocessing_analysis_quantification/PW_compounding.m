%% Cherent and incoherent compounding
%   This example and exercise demonstrates what we in the USTB defines as
%   postprocessing processing objects. Thus, how one can work with the
%   delayed data after beamforming. More specific we are going to look at
%   coherent and incoherent compounding of individual plane wave images.
%   When you do plane wave imaging you can create one low quality image
%   from each transmit. These low quality images can be combined together
%   to create images of higher quality. We will explore this in this
%   exercise.
%
%   The exercise:
%       Part I : Implement both coherent compounding and incoherent
%       compounding of the individual plane wave images. You can read about
%       coherent and incoherent compounding in section 1.7.3 and 1.7.4 in
%       the compendium "Software Beamforming in Medical Ultrasound Imaging"
%
%       Part II : Comparing your implementation to the USTB implementation.
%       
%       Part III : Implement a mix of coherent and incoherent compounding
%
%       Part IV : Compare the resoluting resolution from coherent, 
%                                         incoherent and mix compounding
%
%   Author: Ole Marius Hoel Rindal <olemarius@olemarius.net>

%% Getting the data
%
% We define the local path and the url where the data is stored

% data location
url='http://ustb.no/datasets/';      % if not found data will be downloaded from here
filename='PICMUS_experiment_resolution_distortion.uff';

% checks if the data is in your data path, and downloads it otherwise.
% The defaults data path is under USTB's folder, but you can change this
% by setting an environment variable with setenv(DATA_PATH,'the_path_you_want_to_use');
tools.download(filename, url, data_path);   

%% Loading channel data & scan
%
% The file also contain channel_data and scan. We read it so we can
% replicate the beamformed image in the UFF file.

channel_data=uff.read_object([data_path filesep filename],'/channel_data');
%%
scan=uff.linear_scan()
scan.x_axis = linspace(channel_data.probe.x(1),channel_data.probe.x(end),512)';
scan.z_axis = linspace(0,50e-3,512)';
%uff.read_object([data_path filesep filename],'/scan');

%% Beamforming
%
% We define a pipeline, and the corresponding transmit and apodization
% windows, and launch it.

das = midprocess.das()
das.scan = scan;
das.channel_data = channel_data;
das.dimension = dimension.receive();
das.receive_apodization.window = uff.window.tukey50;
das.receive_apodization.f_number = 1.7;
das.transmit_apodization.window = uff.window.none;
das.transmit_apodization.f_number = 1.7;
b_data_all = das.go()

%% Finally, we will plot the individual plane wave images using the built in
% plot function in USTB
b_data_all.plot([],['Individual PW images'],[],[],[],[],[],'dark')

%% Exercise part I
% Now, in this exercise you are going to explore coherent and incoherent
% compounding of the individual plane wave images. See section 1.7.3 and
% 1.7.4 in Rindal. We wil explore how coherent and incoherent compounding
% of the plane wave images influences the resolution by measuring the width
% of point targets in the image.

% First, let us get all the individual plane wave images in one matrix. The
% "none" argument means that we get the beamformed but complex data. Thus
% the data before envelope detection and log compression.
image_matrix = b_data_all.get_image('none');

% Let us verify that the size of the matrix is 512x512x75 since the number
% of x-pixels is 512, the number of z-axis is 512 and the number of
% transmits, thus individual plane waves are 75
size(image_matrix)

scan.N_x_axis
scan.N_z_axis
channel_data.N_waves

%% Let us define some empty variables to contain the values you calculate
single_image = zeros(scan.N_x_axis,scan.N_z_axis);
coherent_compounding = zeros(scan.N_x_axis,scan.N_z_axis);
incoherent_compounding = zeros(scan.N_x_axis,scan.N_z_axis);
mix_compounding = zeros(scan.N_x_axis,scan.N_z_axis);

% Extract single image
single_image = db(abs(image_matrix(:,:,37)./max(max(image_matrix(:,:,37)))));

% Create coherent compounded image


% Create incoherent compounded image

%% Verify your implementation of coherent and incoherent compounding
% Using changing the dimension to "both" to get coherent compounding
das.dimension = dimension.both()
b_data_coherent = das.go();
b_data_coherent.plot([],['USTB Coherent Compounding using midprocess'])

% Alternatively one can use the coherent compounding postprocess object
cc = postprocess.coherent_compounding();
cc.input = b_data_all;
b_data_coherent = cc.go();
b_data_coherent.plot([],['USTB Coherent Compounding using postprocess'])
coherent_compounding_USTB = b_data_coherent.get_image();

% Using the postprocess incoherent_compunding
ic = postprocess.incoherent_compounding()
ic.input = b_data_all;
b_data_incoherent = ic.go()
b_data_incoherent.plot([],['USTB Incoherent Compounding'])
incoherent_compounding_USTB = b_data_incoherent.get_image();

%% Exercise Part II: Comparing your implementation to the USTB implementation
% Now, you need to plot and verify that your implementation is similar to 
% the images obtained with the USTB. You can do this in a similar matter to
% how you compared your implementation of beamforming in module 3. 
% NB! Do to small numerical differences you can accept a small numerical
% difference between the images and tolerate a pixel difference of
% e.g. abs(img_1_dB - img_2_dB) < 10^-4

%% Compare your implementation of coherent compounding to the USTB

%% Compare your implementation of incoherent compounding to the USTB


%% Exercise Part III : Implement a mix of coherent and incoherent compounding
% We can also to something inbetween full coherent and incoherent
% compounding. We can for example split the low quality images into two,
% and sum the different halfs coherently, and sum those two images
% incoherently. If you for example split the transmit angles into to as:

angles_first_sum = 1:channel_data.N_waves/2;
angles_second_sum = round(channel_data.N_waves/2)+1:channel_data.N_waves;

% And then sum the plane wave images resulting from the angles_first_sum
% coherently, and then angles_first_sum coherently. Before they both are
% combined incoherently.


%% Exercise Part IV : Compare the resoluting resolution from coherent, 
%                                         incoherent and mix compounding
% Below, we provide the code to plot the lateral line, the line along the
% x-axis, through the point scatter at 19 mm. Often resolution is measured
% as the Full Width Half Maximum (FWHM) equal to the width at -6 dB.
% Improved resolution means smaller width of the point scatter. Discuss how
% the different compounding strategies influenced the resolution of this
% point scatter.
line = 191;

figure()
subplot(221)
imagesc(scan.x_axis*1000, scan.z_axis*1000, single_image)
colormap gray; axis image;colorbar; caxis([-60 0])
title('Single transmit');
subplot(222)
imagesc(scan.x_axis*1000, scan.z_axis*1000, coherent_compounding)
colormap gray; axis image;colorbar; caxis([-60 0])
title('Coherent Compounding')
subplot(223)
imagesc(scan.x_axis*1000, scan.z_axis*1000, incoherent_compounding)
colormap gray; axis image;colorbar;  caxis([-60 0])
title('Incoherent Compounding')
subplot(224)
imagesc(scan.x_axis*1000, scan.z_axis*1000, mix_compounding)
colormap gray; axis image;colorbar;  caxis([-60 0])
title('Mix Compounding')


figure();hold all;
plot(scan.x_axis*1000,single_image(line,:),'LineWidth',2,'DisplayName','Single Transmit Image');hold on;
plot(scan.x_axis*1000,coherent_compounding(line,:),'LineWidth',2,'DisplayName','Coherent Compounding');hold on;
plot(scan.x_axis*1000,incoherent_compounding(line,:),'LineWidth',2,'DisplayName','Incoherent Compounding');
plot(scan.x_axis*1000,mix_compounding(line,:),'LineWidth',2,'DisplayName','Mix Compounding');
plot(scan.x_axis*1000,ones(1,scan.N_x_axis)*-6,'r--','LineWidth',2,'DisplayName','- 6dB (FWHM)')
xlim([-3 3]);ylim([-40 0])
legend;
xlabel('x [mm]'); ylabel('Amplitude [dB]')
title(['Lateral line through point scatterer at ',num2str(scan.z_axis(line)*1000,2),' mm']);

