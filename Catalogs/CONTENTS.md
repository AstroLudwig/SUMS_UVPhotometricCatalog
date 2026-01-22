# Contents of the SUMS UV Photometric Catalog
The SUMS UV catalog is split into two tables which each have different magnitude systems:

1. SUMS_UV_Catalog: Provided in _Vega magnitudes_ to be consistent with the original MCPS photometry.
2. SUMS_UV_Catalog_Candidates: Given in _AB magnitudes_.

In all cases the photometry _has not_ been corrected for extinction. The Vega to AB conversions and dust extinction values we use are provided in a table located in the appendix of our paper. <br>

The below table describes columns in the catalog. Columns specific to the _candidate table_ are denoted with an asterisk. 

| Column Name  	| Unit  	| Description  	|   
|---	|---	|---	|
|SUMSID   	| _  	| Unique identifier for each source  	| 
|Galaxy |-|Host galaxy name|
|RAdeg   	|Deg   	|Right Ascension (J2000)   	|   
|DEdeg   	|Deg   	|Declination (J2000)   	|   	
|Ranking* |-|Candidate rank based on degree of UV excess (VB/B) and quality of SED (E/G). [1]|
|UVW2mag |Mag|The UVW2 band magnitude|
|e_UVW2mag |Mag|The UVW2 band magnitude error|
|UVM2mag |Mag|The UVM2 band magnitude|
|e_UVM2mag |Mag|The UVM2 band magnitude error|
|UVW1mag |Mag|The UVW1 band magnitude|
|e_UVW1mag |Mag|The UVW1 band magnitude error|
|Umag |Mag|The MCPS U band magnitude|
|e_Umag |Mag|The MCPS U band magnitude error|
|Bmag |Mag|The MCPS B band magnitude|
|e_Bmag |Mag|The MCPS B band magnitude error|
|Vmag |Mag|The MCPS V band magnitude|
|e_Vmag |Mag|The MCPS V band magnitude error|
|Imag |Mag|The MCPS I band magnitude|
|e_Imag |Mag|The MCPS I band magnitude error|
|UVW2sd |Mag|Standard deviation of the UVW2 band measurements if multiple observations exist|
|UVW2-flux-frac |Mag|Fraction of flux in the UVW2 band within a 5'' radius attributed to the source|
|UVW2-resid-frac |Mag|Fraction of residual flux in the UVW2 band within a 5'' radius after source model subtraction|
|UVW2distm |Arcsec|Distance the source position moved in the UVW2 band during Tractor fitting|
|UVW2distn |Arcsec|Distance to closest neighboring source in the UVW2 band|
|UVW2n5 |Count|Number of sources within 5'' radius in the UVW2 band|
|UVW2n2p5 |Count|Number of sources within 2.5'' radius in the UVW2 band|
|UVW2nobs |Count|Number of observations used in the UVW2 band|
|UVM2sd |Mag|Standard deviation of the UVM2 band measurements if multiple observations exist|
|UVM2-flux-frac |Mag|Fraction of flux in the UVM2 band within a 5'' radius attributed to the source|
|UVM2-resid-frac |Mag|Fraction of residual flux in the UVM2 band within a 5'' radius after source model subtraction|
|UVM2distm |Arcsec|Distance the source position moved in the UVM2 band during Tractor fitting|
|UVM2distn |Arcsec|Distance to closest neighboring source in the UVM2 band|
|UVM2n5 |Count|Number of sources within 5'' radius in the UVM2 band|
|UVM2n2p5 |Count|Number of sources within 2.5'' radius in the UVM2 band|
|UVM2nobs |Count|Number of observations used in the UVM2 band|
|UVW1sd |Mag|Standard deviation of the UVW1 band measurements if multiple observations exist|
|UVW1-flux-frac |Mag|Fraction of flux in the UVW1 band within a 5'' radius attributed to the source|
|UVW1-resid-frac |Mag|Fraction of residual flux in the UVW1 band within a 5'' radius after source model subtraction|
|UVW1distm |Arcsec|Distance the source position moved in the UVW1 band during Tractor fitting|
|UVW1distn |Arcsec|Distance to closest neighboring source in the UVW1 band|
|UVW1n5 |Count|Number of sources within 5'' radius in the UVW1 band|
|UVW1n2p5 |Count|Number of sources within 2.5'' radius in the UVW1 band|
|UVW1nobs |Count|Number of observations used in the UVW1 band|
|Gaia-chi2* |-|Chi-squared statistic from Gaia DR3 proper motions between the source and approximately 1 million members of the host galaxy|

## Notes:
[1] Candidate rankings are a combination of UV excess (VB/B) and SED quality (E/G) flags, described below: <br>
- Very Blue (VB): More than 0.4 mag bluer than the ZAMS at the source's brightness. <br>
- Blue (B): Up to 0.4 mag bluer than the ZAMS at the source's brightness. <br>
- Excellent (E): Average UV flux fraction is >40%, at least 6 of the 7 photometric bands present, SED is highly consistent with models. <br>
- Good (G): Average UV flux fraction is >10% if closest source is further than 2.5'' otherwise > 25%, at least 5 of the 7 photometric bands present, SED is consistent with models.  <br>

  _Additional details on these rankings can be found in Section 5 of our paper._ <br>
