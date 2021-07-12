# Module 7 : Advanced Methods in Ultrasound Imaging

In this exercise you will exlore more advanced beamforming, often known
as adaptive beamforming. Because we are adapting the beamforming process
to the data we have received.

More specifically, we will look at the coherence factor.

## Part I:
    Calculate the coherence factor as in equation 1.38 in the compendium.

## Part II:
    Check your implementation agains the USTB implementation.

## Part III:
    Analyse the delayed data:
    + How will your observations of delay affect the coherence factor?
    + What is the image you plot of the coherent sum equal to?

## Part IV:
    Applying the CF as a image weight to the DAS image

## Part V:
    Compare DAS CF to DAS image:
    + What are the differences between the results of the conventional DAS to the image with DAS weighted with CF?
    + What happened to the object from x = 0 to x = 2.5 mm at z = 30 mm in the two plots?
    + In the plot below we also plot the mean lateral line through the gradient
        from x = +-14mm at z = 40 to 48 mm. Theoretically, this should go from 0
        to -50 dB, which one is most correct?
