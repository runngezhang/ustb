# Module 2 : Wave Physics

This example demonstrates how to run a simulation in k-wave recording a
single originating from a single source. We are then showing how we can
import the recorded signal into a channel data UFF object in the USTB and
beamforming the data into an image. Thus we will do "receive beamforming"
and reconstruct an image of the single source. Your task is to implement
your own pixelbased receive beamformer.

## Litterature:
For background you can read page 1 and 2 of Jørgen Grythes document
"Beamforming Algorithms - beamformers" or pages 22-29. However, remember
that here you only need to do receive beamforming. The compendium for the
course is relevant, especially section 1.7.

## The exercise:
### Part I
Implement your own pixelbased beamformer.

### Part II
Reflect and answer the following questions:
+ What happens when you change the number of sensors from 4 to 16?
+ What happens when you change the transmit signal from *gausian_pulse* to *sinus*?
+ What is illustrated in Figure 10? Explain the images and how they differ from the final image.

### Part III
Visualize the channel data before and after delay for point scatter
First of all, this plot is much better if you use e.g. 16 elements on
line 32. So you should go ahead and change that. 

Your task here is to use the plot above to find the location of the point
scatter. Use the cursor in the plot and find the maximum, and simply set
the correct value in the variables below for the x and z locatino of the
scatterer, also known as the point spread function (PSF).
