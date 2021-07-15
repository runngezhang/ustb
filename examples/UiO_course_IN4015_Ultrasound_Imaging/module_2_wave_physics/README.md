# Module 2 : Wave Physics

This example demonstrates how to run a simulation in k-wave recording a
single originating from a single source. We are then showing how we can
import the recorded signal into a channel data UFF object in the USTB and
beamforming the data into an image. Thus we will do "receive beamforming"
and reconstruct an image of the single source. Your task is to implement
your own pixelbased receive beamformer.

## Litterature:
Page 1 and 2 of Jørgen Grytes "Beamforming Algorithms - beamformers" or
pages 22-29 in the compendium. 

## The exercise:
### Part I
Implement your own pixelbased beamformer.

### Part II
Reflect and answer the following questions:
+ What happens when you change the number of sensors from 4 to 16?
+ What happens when you change the transmit signal from *gauusian_pulse* to *sinus*?
+ What is illustrated in Figure 10? Explain the images and how they differ from the final image.
