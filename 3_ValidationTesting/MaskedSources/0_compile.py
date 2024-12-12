import glob 
import pandas as pd 
import numpy as np
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u

# First compile catalog of all masked sources 
# Get gaia / simbad crossmatch in topcat 
# Remove duplicates 

for galaxy in ['lmc','smc']:

    # FILE HANDLING
    files = glob.glob(f"H:/Data/SUMS_Tractor_Data/{galaxy}/*X/*/*masked*")

    # Combine all files 
    df = pd.concat([pd.read_csv(f) for f in files])
    init_size = df.shape[0]

    # Sort by RA
    df = df.sort_values('ra')
    df = df.reset_index(drop=True)

    # Find duplicate sources if ra/dec is equivalent, will do additional step in topcat 
    df = df.drop_duplicates(subset=['ra','dec']).reset_index(drop=True)

    df.to_csv(f'C:/Projects/0_Data/SUMS_CompleteCatalog/MaskedSources/{galaxy}_masked_sources.csv',index=False)