# Caveats and Limitations of the SUMS UV Photometric Catalog
There are a number of caveats and limitations to be aware of when using this catalog. 
Here, we provide a short summary for a number of these issues, as well as flags that are included in the final catalog to help assess and mitigate their effects.
Further details can be found in the Appendix of our paper and the sections referenced therein. 

## :pushpin: Incomplete Source Coverage:
This UV catalog is not complete, even for sources above the detection threshold. 

### Extremely Bright Sources and their Surroundings are Masked
We removed sources within 12" of stars bright enough to cause aberrations that would corrupt photometric estimates in their local environment. 
This bright source masking impacts <1% of the total coverage for each galaxy. 

### Densely Clustered Regions Have Higher Photometric Uncertainty
We are likely less complete in densely clustered regions, and regions close to bright (but un-masked) stars. 
We remove sources with high residual flux fractions to account for associated photometric errors, but this means we are likely less sensitive to faint sources in these regions.

### Unmodeled Sources:
Only sources in the MCPS catalog were modeled by the Tractor. The drift scan used in MCPS includes gaps which therefore propogate to our catalog. 
Further, for a source from MCPS to be included in the input list for the Tractor to model, we require that it is brighter that 20.5 mag in either the U- or B-band 
and that aperture photometry within a 5" radius yields a count rate that is >1.5 times higher than the background model in the same region. 
This means that some sources for which there is flux in the UVOT images do not appear in the final catalog because they were not modeled. 

_This repo includes our photometric pipeline for any user who would like to examine the outputs for particular targets of interest that are not within our catalog._

## :pushpin: Flux Degeneracy for Close Sources:
Degeneracies can occur between the count rates that the Tractor assigns to sources that are physically close together. 

In our mock-source injection test, the divergence begins to appear once stars are closer than ~2-2.5'' (i.e. within the FWHM of the Swift-UVOT PSF and is more extreme for the fainter object. 

_To help identify nearby source flux misattribution, the catalog includes the fraction of the flux within a 5" radius attributed to a given source, the distance to the nearest neighbor, and the number of neighbors within 2.5 and 5"._

## :pushpin: PSF Shape Variations

We use a constant theoretical model of the Swift-UVOT PSF, however, there are two regimes that may be negatively impacted by this: 

- For bright stars, coincidence narrows the PSF which could impact the magnitudes we get for those stars. 
- Although we remove images with clear tracking issues, some images may still vary from the ideal case represented by the model PSF. 

To account for this, we add a 5% systematic error to the measured count rates computed with the Tractor and we validate our pipeline against standard Swift routines in less crowded regions. 

_We also establish the `residual fraction' column as an indicator of how the flux of a source could have been underestimated, therefore suggesting that the fit may not have been as robust._ 

## :pushpin: UV-Optical Source Mismatches
We calculate UV magnitudes by performing forced-photometry at the location of MCPS sources. However, the MCPS source and the Swift-UVOT source may become mismatched for the following reasons: 

- The MCPS catalog is not complete at bright magnitudes, which may result in a bright optical source being matched to a nearby source in the UV image. This would become apparent in the SED when a large discrepancy between the UV and optical magnitudes occurs.  
- The Tractor pipeline fits the position of a source up to 1'' of the initial position to account for astrometric errors. In densely clustered regions, this could result in a mismatch which would similarly result in an SED discontinuity.

_In addition to the nearest neighbor columns, we include the 'dist_moved' columns as a way of determining how far the Tractor optimized the position from the initial MCPS coordinate._ 

