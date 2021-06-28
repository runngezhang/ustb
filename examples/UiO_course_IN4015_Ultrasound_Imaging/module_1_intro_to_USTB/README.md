# Module 1 Introduction to the USTB: Recording and prosessing ultrasound images

The lab assignments for module 1 consists of using existing USTB examples for part 2 and 3 of the exercise.

+ Part 2: Get to know the UltraSound ToolBox by running and getting familiar with multiple examples 
	+ Pure USTB examples. Run at least three of these
		+ examples/uff/CPWC_UFF_Alpinion.m
		+ examples/uff/CPWC_UFF_Verasonics.m
		+ examples/uff/FI_UFF_phased_array.m
		+ examples/picmus/experimental_contrast_speckle.m
		+ examples/picmus/experimental_resolution_distortion.m
		+ examples/picmus/carotid_cross.m
		+ examples/picmus/carotid_long.m
		+ examples/acoustical_radiation_force_imaging/ARFI_UFF_Verasonics.m
        + examples/UiO_course_IN4015_Ultrasound_Imaging/module_1_intro_to_USTB/minimal_example.m
        + examples/UiO_course_IN4015_Ultrasound_Imaging/module_1_intro_to_USTB/maximal_example.m
	+ USTB + K-wave examples. You need to install and add the ultrasound simulator k-wave (http://www.k-wave.org/) to your MATLAB path. See simple description below.
		+ examples/kWave/CPWC_linear_array_cyst.m
	+ USTB + Field II examples. You need to install and add ultrasound simulator Field II (https://field-ii.dk/) to your MATLAB path. See simple description below.
		+ examples/field_II/STAI_L11_resolution_phantom.m
		
+ Part 3 : Recording an image of a CIRS phantom using the Verasonics Vantage 256 scanner. 
	+ Using the example available at examples/verasonics/FI_linear_array_L11.m record two images
		+ One image of the point scatteres in the phantom
		+ One image of the cysts in the phantom
	+ Using the example available at examples/verasonics/FI_phased_array_P4.m record an image covering both the point scatteres and the cysts in the phantom.

## How to download k-wave and Field II

### K-wave

Register a new profile (http://www.k-wave.org/forum/register.php) with the necessary information, go to download, and download the k-wave zip-file.
Next, extract the folder on your preferred spot, and add to path (see below).

### Field II

Go to the download site at (https://field-ii.dk/), and download the file for your operating system. Extract the folder (twice) at your preferred location.
Next, add to path (see below).

### Add to path

There are two ways to add to path:

1. Create a new file called "startup.m", and add the path through this file by writing
    + addpath *filepath, i.e. Documents/MATLAB/k-wave*
2. Pathtool: open pathtool in Matlab by writing pathtool in the Command Window in Matlab. Click "Add Folder" to add k-wave and field-ii to path, and press save. To use this method, you have to have Admin rights on the computer you are using.
