# Contents of the SUMS UV Photometric Catalog
The SUMS UV catalog is split into three tables:

1. SUMS_UV_Catalog_LMC
2. SUMS_UV_Catalog_SMC
3. SUMS_UV_Catalog_Candidates

The first two tables are provided in _Vega magnitudes_ to be consistent with the original MCPS photometry. <br>
However, the candidate catalog (table 3) is given in _AB magnitudes_. <br>

In all cases the photometry _has not_ been corrected for extinction. The Vega to AB conversions and dust extinction values we use are provided in Table A1 of our paper. <br>

The below table describes columns in the catalog. Columns specific to the _candidate table_ are denoted with an asterisk. 

| Column Name  	| Unit  	| Description  	|   
|---	|---	|---	|
|SUMS_ID   	| _  	| Unique identifier for each source  	| 
|Galaxy* |-|Host galaxy name|
|RA   	|Deg   	|Right Ascension (J2000)   	|   
|Dec   	|Deg   	|Declination (J2000)   	|   	
|Ranking* |-|Candidate rank based on degree of UV excess (VB/B) and quality of SED (E/G). [1]|
|UVW2 |Mag|The UVW2 band magnitude|
|UVW2_err |Mag|The UVW2 band magnitude error|
|UVM2 |Mag|The UVM2 band magnitude|
|UVM2_err |Mag|The UVM2 band magnitude error|
|UVW1 |Mag|The UVW1 band magnitude|
|UVW1_err |Mag|The UVW1 band magnitude error|
|U |Mag|The MCPS U band magnitude|
|U_err |Mag|The MCPS U band magnitude error|
|B |Mag|The MCPS B band magnitude|
|B_err |Mag|The MCPS B band magnitude error|
|V |Mag|The MCPS V band magnitude|
|V_err |Mag|The MCPS V band magnitude error|
|I |Mag|The MCPS I band magnitude|
|I_err |Mag|The MCPS I band magnitude error|
|UVW2_std |Mag|Standard deviation of the UVW2 band measurements if multiple observations exist|
|UVW2_flux_frac |Mag|Fraction of flux in the UVW2 band within a 5'' radius attributed to the source|
|UVW2_resid_frac |Mag|Fraction of residual flux in the UVW2 band within a 5'' radius after source model subtraction|
|UVW2_dist_moved |Arcsec|Distance the source position moved in the UVW2 band during Tractor fitting|
|UVW2_dist_neighbor |Arcsec|Distance to closest neighboring source in the UVW2 band|
|UVW2_n5 |Count|Number of sources within 5'' radius in the UVW2 band|
|UVW2_n2p5 |Count|Number of sources within 2.5'' radius in the UVW2 band|
|UVW2_nobs |Count|Number of observations used in the UVW2 band|
|UVM2_std |Mag|Standard deviation of the UVM2 band measurements if multiple observations exist|
|UVM2_flux_frac |Mag|Fraction of flux in the UVM2 band within a 5'' radius attributed to the source|
|UVM2_resid_frac |Mag|Fraction of residual flux in the UVM2 band within a 5'' radius after source model subtraction|
|UVM2_dist_moved |Arcsec|Distance the source position moved in the UVM2 band during Tractor fitting|
|UVM2_dist_neighbor |Arcsec|Distance to closest neighboring source in the UVM2 band|
|UVM2_n5 |Count|Number of sources within 5'' radius in the UVM2 band|
|UVM2_n2p5 |Count|Number of sources within 2.5'' radius in the UVM2 band|
|UVM2_nobs |Count|Number of observations used in the UVM2 band|
|UVW1_std |Mag|Standard deviation of the UVW1 band measurements if multiple observations exist|
|UVW1_flux_frac |Mag|Fraction of flux in the UVW1 band within a 5'' radius attributed to the source|
|UVW1_resid_frac |Mag|Fraction of residual flux in the UVW1 band within a 5'' radius after source model subtraction|
|UVW1_dist_moved |Arcsec|Distance the source position moved in the UVW1 band during Tractor fitting|
|UVW1_dist_neighbor |Arcsec|Distance to closest neighboring source in the UVW1 band|
|UVW1_n5 |Count|Number of sources within 5'' radius in the UVW1 band|
|UVW1_n2p5 |Count|Number of sources within 2.5'' radius in the UVW1 band|
|UVW1_nobs |Count|Number of observations used in the UVW1 band|
|Gaia_chi2* |-|Chi-squared statistic from Gaia DR3 proper motions between the source and approximately 1 million members of the host galaxy|

## Notes:
[1] Candidate rankings are a combination of UV excess (VB/B) and SED quality (E/G) flags, described below: <br>
- Very Blue (VB): More than 0.4 mag bluer than the ZAMS at the source's brightness. <br>
- Blue (B): Up to 0.4 mag bluer than the ZAMS at the source's brightness. <br>
- Excellent (E): Average UV flux fraction is >40%, at least 6 of the 7 photometric bands present, SED is highly consistent with models. <br>
- Good (G): Average UV flux fraction is >10% if closest source is further than 2.5'' otherwise > 25%, at least 5 of the 7 photometric bands present, SED is consistent with models.  <br>

  _Additional details on these rankings can be found in Section 5 of our paper._ <br>
