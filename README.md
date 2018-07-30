# EcoFOCI_Prawler
Collection of tools to work with PMEL-Engr Prawler data from the yawl/ketch data processing system (this data has had engineering calibrations applied and is a clean organized human readable version of the rudics data)

### Gridding Procedure

Data is distributed as ragged arrays where each profile may have slightly varying lengths based on number of samples and the fall speed.  Adding the park and hold calibration samples compounds this issue as each hold is considered a _dive_ .

If gridded data is desired a few options need to be addressed.  Below are the basic steps/assumptions.
**Profile**
    + Grid in the vertical: values greater than 0.25m seem acceptable.  
    + If multiple data points exist in this depth bin, take the median value. (`press_grid`)
    + Linearly fill in between depth bins when no data is present (`fill gaps = True`)
 
 *
 **In Time**
    + We're already gridded in depth so set up a grid in time (hourly is appropriate)
    + Use Python-Scipy 2d interpolate function which is a linear spline interpolation (B-spline of a surface) and set fill values to -100000 if outside domain range


################

Legal Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce (DOC). All NOAA GitHub project code is provided on an 'as is' basis and the user assumes responsibility for its use. Any claims against the DOC or DOC bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation, or favoring by the DOC. The DOC seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by the DOC or the United States Government.