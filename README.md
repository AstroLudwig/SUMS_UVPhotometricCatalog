# The Stripped-Star Ultraviolet Magellanic Cloud Survey (SUMS):
The UV Photometric Catalog and Stripped Star Candidate Selection

![The Complete Catalog](/2_CompileCatalog/sums_photometry.png "SUMS Catalog")
# Authors & License
Copyright XX-XX Bethany Ludwig (KU Leuven), Maria Drout (University of Toronto), Ylva Gotberg (Institute of Science and Technology Austria), Dustin Lang (Perimeter Institute), and Alexander Laroche (University of Toronto).<br>

Licensed under XXXX; see [LICENSE](XXX).<br>

Please contact the first author should you have any questions.<br>

# Abstract & Summary
This repository contains all code used in this project. Below is a brief summary, please consult our paper for further details.<br> 

Most massive stars (~8-25 Msun) interact with a binary companion during their lifetimes. These interactions can remove the hydrogen-rich envelope, producing intermediate-mass (~2-8 Msun) and helium-rich stars. These **stripped stars** are predicted to emit predominantly in the ultraviolet (UV) and can therefore be identified via a UV excess, provided they are not outshone by their companion. However, despite their importance to binary evolution, supernovae, and ionizing feedback, few stripped stars have been confirmed. This is likely due to the scarcity of wide-field, high angular-resolution, UV surveys of stellar populations with reliable distances and extinction estimates. To address this, we present the Stripped-Star Ultraviolet Magellanic Clouds Survey (SUMS) catalog.<br>

In 0_CleanData we obtain archival Swift-UVOT images of the LMC and SMC. [Astrometry](https://astrometry.net/) with [custom index files](https://zenodo.org/records/7600727/files/S311_astrometry_index_files.zip?download=1) is run on these files to update their WCS.<br>
   
In 1_TractorPipeline we use [the Tractor](https://github.com/dstndstn/tractor) forward modeling software to perform PSF photometry on these images. The photometry is converted to magnitudes using the standard Swift-UVOT calibration tools provided by [HEASARC](https://heasarc.gsfc.nasa.gov/docs/software.html) with the following calibration files:<br>
- [swusenscorr20041120v005.fits](https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/swift/docs/uvot/uvot_release_history.html)
- [sssfile5.fits.gz](https://swift.gsfc.nasa.gov/analysis/uvot_digest/sss_check.html)
  
In 2_CompileCatalog we assemble the catalog, resulting in 734,862 sources in three UV filters to a depth of ~20 Vega mag. The full SUMS UV catalog for the LMC and SMC can be found in the Catalogs folder. Please note that these magnitudes are in Vega with no extinction correction. A table with the AB conversion and our assumed extinction correction can be found in 0_CleanData and in the appendix of our paper. <br>
  
In 3_ValidationTesting we perform a series of validation tests on the photometry pipeline.<br>  
  
In 4_FindCandidates we identify sources with excess UV light compared to main-sequence stars and apply a series of quality and SED cuts. From this, we identify 522 candidate stripped stars in the LMC and 298 in the SMC. The candidate catalog can also be found in the Catalogs folder, however, note that these magnitudes are in AB.<br>  
  
In 5_AssessCandidates we assess the properties of our candidates including their overall brightnesses, SEDS, and spatial distribution.<br>   
   
In 6_OtherBlue we assess the potential contamination from early main-sequence stars and other UV excess systems including foreground and background sources.<br>  
   
This survey lays the groundwork for the first systematic census of stripped stars and opens new windows into binary evolution and massive star populations.  


