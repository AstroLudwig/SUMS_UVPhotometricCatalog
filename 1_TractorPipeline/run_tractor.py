import sys
# Gave Ylva code to run on cluster since it takes a long time. 
sys.path.insert(0, '/data002/ygoetberg/Tractor/astrometry.net-0.85/')

import matplotlib.pyplot as plt 
import numpy as np 
import os,time
import glob,sys

from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable

sys.path.insert(0, 'lib/')
from RetrieveSource import *
from EstimateBackground import *
from TractorTools import *
from PSF import *

# Input ########################################################################################################
# UVOT Filters 
filters = ["um2","uw1","uw2"]

# Observations Ids for SMC_# Survey
smc_observation_id = np.arange(40415,40465,1).astype(str)

# Observations Ids for LMC_Survey_# Survey
lmc_observation_id = np.arange(45422,45587,1).astype(str)

uvfilter = filters[2]   
obsid = 45490
galaxy = 'lmc'
segment = 1
extension = 2 

path = 'data'
###########################################################################################################################

# Get Files
# This is the MPCS catalog for this field. Always ends in .full.dat
cname = f'{path}/sw000{obsid}00{segment}{uvfilter}_sk_{obsid}_{segment}_{extension}.full.dat' 
# This is the data. Always ends in .new
fname = f'{path}/sw000{obsid}00{segment}{uvfilter}_sk_{obsid}_{segment}_{extension}.new' 
# Where to save
directory = f"photometry/{galaxy}/{obsid}/"
# What to save output as. All the other files that are saved are based on this name.
savefile = directory + f"{obsid}_{uvfilter}_{segment}_{extension}.fits" 

# Create directory
if not os.path.exists(directory):
    os.makedirs(directory)

# Get Catalog
labels = ['RAhr','DEdeg','Umag','e_Umag','Bmag','e_Bmag','Vmag','e_Vmag','Imag','e_Imag','Flag','Jmag','e_Jmag','Hmag','e_Hmag','Ksmag','e_Ksmag']
cat = pd.read_csv(cname,delimiter='\s+',names=labels)

# Get Data
hdr = fits.open(fname)[0]

# Timer for reference
start = time.time()

print('---Getting Catalog---')
source = get_meta().with_hdu(hdu=hdr,
                 usno_catalog=f'usno/{galaxy}_anti_match_USNO_inf.dat', # For removing bright things not in mcps
                 optical_catalog=cname,
                 directory=directory,
                 Umag_cutoff=20.5,
                 Bmag_cutoff=20.5,
                 fits_origin=0,
                 aperture_size=2.5*2,
                 xdim=[0,np.shape(hdr.data)[1]],  #  <---- Can change how much of the field you want to run on. Makes runs much shorter.
                 ydim=[0,np.shape(hdr.data)[0]],
                 #xdim=[500,700],  #  <---- Reasonable numbers for quick test.
                 #ydim=[500,700],            
                 save_dropped_catalog=True)


print('---Getting Background---')
bkgd = BkgdEstimator(source,n_pix = [20,20])


print('---Getting PSF---')
pix_scale = np.abs(source.cdelt)*3600.

psf_object = psf_fit(pixel_per_arsecond = 1/pix_scale,
                     uvfilter = uvfilter, width = 23,
                     cog_file='calib/swureef20041120v104.fits').psf


print('---Running Tractor---')

TractorObject = PhotometryTools(source,
                    psf_filename = psf_object,
                    fits_table_sname = savefile, 
                    background = bkgd,
                    fit_positions = np.nan, # This is what takes so long, you can change it to np.nan for a quick test. Otherwise it's 0.05
                    threshold = 1.5,
                    save_output=True)


print(f"Time Taken for Tractor: {(time.time()-start)/60/60} hours")

