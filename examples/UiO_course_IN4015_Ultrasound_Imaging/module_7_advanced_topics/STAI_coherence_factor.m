% % Exercise for Module 7 : Advanced Methods in Ultrasound Imaging
% In this exercise you will exlore more advanced beamforming, often known
% as adaptive beamforming. Because we are adapting the beamforming process
% to the data we have received.
% 
% More specifically, we will look at the coherence factor.
% See the parts below.
% 
% Part I:
%     Calculate the coherence factor as in equation 1.38 in the compendium
% Part II:
%     Check your implementation agains the USTB implementation
% Part III:
%     Analyse the delayed data
% Part IV:
%     Applying the CF as a image weight to the DAS image
% Part V:
%     Compare DAS CF to DAS image
%
% The dataset might fail to download, if so delete it from the folder data,
% and download it from this Google drive link:
% https://drive.google.com/file/d/1xAXoEWhPcYjam9R1iuQ0gWKdDXiVlPCX/view?usp=sharing


filename = 'FieldII_STAI_dynamic_range.uff';
url='http://ustb.no/datasets/';      % if not found downloaded from here

% Download data
tools.download(filename, url, data_path);   

% Read data
channel_data = uff.channel_data();
channel_data.read([data_path,filesep,filename],'/channel_data');

%% Define scan
scan=uff.linear_scan('x_axis',linspace(-20e-3,20e-3,512).', 'z_axis', linspace(10e-3,52.5e-3,512).');

%% Beamformer
mid = midprocess.das();
mid.channel_data=channel_data;
mid.scan=scan;
mid.dimension = dimension.transmit();

mid.receive_apodization.window=uff.window.boxcar;
mid.receive_apodization.f_number=1.75;

mid.transmit_apodization.window=uff.window.boxcar;
mid.transmit_apodization.f_number=1.75;
b_data_tx = mid.go();

%% Conventional DAS image
mid.dimension = dimension.both()
b_data_das = mid.go();

%%
b_data_das.plot([],[],[80])
%%
data_cube = reshape(b_data_tx.data,scan.N_z_axis,scan.N_x_axis,channel_data.probe.N_elements);

% To get the scalar M in equation 1.38 in the compendium you can do this.
% M is thus dependent on the receive apodization because different number
% of receive elements are used for different depths (expanding aperture)
w_rx = reshape(mid.receive_apodization.data,scan.N_z_axis,scan.N_x_axis,channel_data.probe.N_elements);
M = sum(w_rx,3);

%% Let us define some empty variables to contain the values you calculate
coherent_sum = zeros(scan.N_x_axis,scan.N_z_axis);
incoherent_sum = zeros(scan.N_x_axis,scan.N_z_axis);
CF = zeros(scan.N_x_axis,scan.N_z_axis);
das_weighted_CF = zeros(scan.N_x_axis,scan.N_z_axis);

%% Exercise Part I: Calculate the coherence factor as in equation 1.38 in the compendium
% Calculate the coherent sum over the elements (the expression in the
% numerator (above the brøkstrek ;))
coherent_sum;
% Calculate the incoherent sum over the elements (the sum in the expression 
% in the denominator (below the brøkstrek ;))
incoherent_sum;
% Use the coherent, incoherent sum and the scalar M for each pixel to
% calculate the coherence factor as in expression 1.38.
CF;

%% USTB implementation of coherence factor
cf = postprocess.coherence_factor()
cf.transmit_apodization = mid.transmit_apodization;
cf.receive_apodization = mid.receive_apodization;
cf.input = b_data_tx;
b_data_cf = cf.go()
das_CF_USTB = b_data_cf.get_image();
USTB_coherence_factor = cf.CF.get_image('none');

%% Part II : Check your implementation against the USTB implementation


%% Part III : Let's analyse the delayed data
% If we plot the delayed data for two pixels in the image, and deliberately
% choosing to plot the delayed data around the point scatter at x = -5.5 mm 
% and z = 35 +- 2 mm, the top plot, and right next to it in x = -4.7 mm with 
% the same z = 35 +- 2mm in the plot below.

figure
subplot(211)
imagesc(18:77,scan.z_axis(302-20:302+20)*1000,squeeze(real(data_cube(302-20:302+20,186,18:77))))
line([18,77],[scan.z_axis(302)*1000,scan.z_axis(302)*1000],'color','r','LineWidth',2)
colorbar
xlabel('Rx Element')
ylabel('z [mm]')
title(['The delayed data through x= ',num2str(scan.x_axis(186)*1000,2),'mm , z = ',num2str(scan.z_axis(302)*1000,2),' mm']);
subplot(212)
imagesc(18:77,scan.z_axis(302-20:302+20)*1000,squeeze(real(data_cube(302-20:302+20,197,18:77))))
line([18,77],[scan.z_axis(302)*1000,scan.z_axis(302)*1000],'color','r','LineWidth',2)
colorbar
xlabel('Rx Element')
ylabel('z [mm]')
title(['The delayed data through x= ',num2str(scan.x_axis(197)*1000,2),'mm , z = ',num2str(scan.z_axis(302)*1000,2),' mm']);

% In the top plot we see that the delayed data align perfectly along the 
% red line overlayed at the location of the point scatter. While in the
% plot below, we see that the the data does not align perfectly and that
% the amplitude of the data is alternating from a positive to a negative
% signal value along the overlaid line.
%
% How will these observations affect the coherence factor?

%% Plotting the coherent and incoherent sum as independent images
figure
subplot(121)
imagesc(scan.x_axis*1000, scan.z_axis*1000, db(abs(coherent_sum./max(coherent_sum(:)))));
axis image; title('Image of coherent sum'); xlabel(['x [mm]']); ylabel(['y [mm]']);
caxis([-200 0])
subplot(122)
imagesc(scan.x_axis*1000, scan.z_axis*1000, db(abs(incoherent_sum./max(incoherent_sum(:)))));
axis image; title('Image of incoherent sum'); xlabel(['x [mm]']); ylabel(['y [mm]']);
colormap gray
caxis([-200 0])

% What is the image of the coherent sum equal to?

%% Part IV: Applying the CF as a image weight to the DAS image
% The default value from the postprocess Coherence Factor was the
% conventional DAS image multiplied with the CF as an image weight as in
% equation 1.39. Implement this yourself and verify that it is equal to the
% USTB implementation.
[weights,array_gain_compensation,geo_spreading_compensation] = tools.uniform_fov_weighting(mid);
das_img_signal = b_data_das.get_image('none').*weights;
das_img_db = db(abs((das_img_signal)./max(das_img_signal(:))));

das_weighted_CF = das_img_signal.*CF;
das_weighted_CF = db(abs(das_weighted_CF./max(das_weighted_CF(:))));

% Compare that to the USTB versions as well. It is in das_CF_USTB


%% Part V: Compare DAS CF to DAS image
% Now, let's compare the results from conventional DAS to the image with
% DAS weighted with CF. What are the differences?
% What happened to the object from x = 0 to x = 2.5 mm at z = 30 mm?
% In the plot below we also plot the mean lateral line through the gradient
% from x = +-14mm at z = 40 to 48 mm. Theoretically, this should go from 0
% to -50 dB, which one is most correct?
%
% If you want to understand why this is happening, you have to take our
% next course! Or read the reference
channel_data.print_authorship

%%
[~,x_start] = min(abs(scan.x_axis+14/1000));
[~,x_stop] = min(abs(scan.x_axis-14/1000));
%%
figure
subplot(221)
imagesc(scan.x_axis*1000, scan.z_axis*1000, das_img_db)
colorbar; caxis([-60 0]); colormap gray; axis image; title('DAS');
subplot(222)
imagesc(scan.x_axis*1000, scan.z_axis*1000, das_weighted_CF)
colorbar; caxis([-60 0]); colormap gray; axis image; title('DAS weighted with CF');
subplot(2,2,[3 4]);hold on;
plot(scan.x_axis(x_start:x_stop)*1000,linspace(0,-50,x_stop-x_start+1),'k--','LineWidth',2,'Displayname','DAS')
plot(scan.x_axis(x_start:x_stop)*1000,mean(das_img_db(350:470,x_start:x_stop),1)-max(mean(das_img_db(350:470,x_start:x_stop),1)),'b','LineWidth',2,'Displayname','DAS weighted CF')
plot(scan.x_axis(x_start:x_stop)*1000,mean(das_weighted_CF(350:470,x_start:x_stop),1)-max(mean(das_weighted_CF(350:470,x_start:x_stop),1)),'r','LineWidth',2,'Displayname','Theoretical')
xlabel('x [mm]');ylabel('Amplitude [dB]');legend show