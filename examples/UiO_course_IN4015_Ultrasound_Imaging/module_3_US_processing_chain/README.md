# Module 3 : US Processing Chain

The lab assignment for **Module 3** is to implement your own beamforming from
scratch. In Module 2, we used pixel based beamforming, calculating the
propagation delay one way, the recieve delay. In this module, we are taking
both the receive delay and the transmit deay of the total propagation delay
into account, using beambased beamforming. This means that we are now using
angles and distance to find the positions, rather than the predefined pixels
as we did in Module 2.

After your implementation of the total propagation delay, your method will
be compared to USTB built in method. In this module we are aiming to minimize
the difference between your method and USTBs method. Finally we will implement
your method to make your own image.


## Litterature
The relevant litterature for this assignment is the lecture slide. Especially the
slides on beambased beamforming. However, the compendium: “Software Beamforming
in Medical Ultrasound Imaging”, is relevant. Especially section 1.5 to 1.7 as 
well as 1.9. However, the other sections prior to 1.5 is also very relevant.
But remember that you are implementing a beambased beamformer, while the compendium
describes the pixel based beamformer in the USTB.

%% Datasets
You have two awailable datasets you can use for this exercise

+ Verasonics_P2-4_parasternal_long_small.uff which is a in-vivo cardiac dataset
+ FieldII_P4_point_scatterers.uff which is a simulated dataset of point scatters

It is perhaps easiest to get your code running correctly using the dataset
with the point scatterers. However, you should also try to use the cardiac dataset.

+ You will also have the option to record your own dataset on the Verasonics scanner
in the lab. More information regardin this will be given in the lecture and the group lecture


## The exercise:
### Part I : Do phased array beamforming with the USTB

Here you dont't have to implemen anything. Just run the code as it is.
You will later compare your beamformed image with the image resulting
from this beamforming. Just try to understand what is going on.

### Part II : Implement your own beambased beamforming from scratch.

Please see the slides on beambased beamforming, and especially the slide
on the geometry og beambased beamforming on tips on how to implement this.

A hint is to review the exercise from module 2 on wave physics where you
implemented an receive beamformer. Now you are extending it to also
include to compensate for the transmit part of the propagation delay as
well as handling multiple transmits.

Pseudocode for beambased phased array beamformer
for each transmit
    calculate the transmit part and the receive part of the delay for
    every receive channel. Remember to calculate the delays in seconds
    not distance, and remember to subtract the offset for each transmit
    event to get a correct time zero convention.

    for each receive channel
        beamform by interpolating (using interp1) into the calculated
        delays for each radial dept sample for each transmit
        the call to interp1 might look like this:
        interp1(sample_time, rfData(:,r,t)', delays(t,:,r))'
            where r is the current receive channel and t is the current transmit
            you allready have sample_time and rfData, but need to calculate
            the propriate delays

