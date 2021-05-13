# Module 1 Introduction to the USTB: Recording and prosessing an ultrasound image

The lab assignments for module 1 consists of using existing USTB examples for part 2 and 3 of the exercise.

+ Part 2 : Recording an image of a CIRS phantom using the Verasonics Vantage 256 scanner. 
	+ Using the example available at examples/verasonics/FI_linear_array_L11.m record two images
		+ One image of the point scatteres in the phantom
		+ One image of the cysts in the phantom
	+ Using the example available at examples/verasonics/FI_phased_array_P4.m record an image covering both the point scatteres and the cysts in the phantom.
	
+ Part 3: Get to know the UltraSound ToolBox by running and getting familiar with multiple examples 
	+ Pure USTB examples. Run at least three of these
		+ examples/picmus/carotid_cross.m
		+ examples/picmus/carotid_long.m
		+ examples/picmus/experimental_contrast_speckle.m
		+ examples/picmus/experimental_resolution_distortion.m
		+ examples/uff/CPWC_UFF_Alpinion.m
		+ examples/uff/CPWC_UFF_Verasonics.m
		+ examples/uff/FI_UFF_phased_array.m
		+ examples/acoustical_radiation_force_imaging/ARFI_UFF_Verasonics.m
	+ USTB + K-wave examples. You need to install and add the ultrasound simulator k-wave (http://www.k-wave.org/) to your MATLAB path.
		+ examples/kWave/CPWC_linear_array_cyst.m
	+ USTB + Field II examples. You need to install and add ultrasound simulator Field II (https://field-ii.dk/) to your MATLAB path.
		+ examples/field_II/STAI_L11_resolution_phantom.m
		
		