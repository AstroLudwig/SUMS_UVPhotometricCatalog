import glob 
import os
import pandas as pd 
import numpy as np
from astropy.coordinates import SkyCoord

test1_step1files_completed = False
test2_step3_uniquecoords_areunique = False
test3_step3_uniquecoords_assignment = True

# Files
save_dir = "C:/Projects/0_Data/SUMS_CompleteCatalog/"
step2_dir = "C:/Projects/0_Data/SUMS_CompleteCatalog/Step2/"
log_dir = "C:/Projects/0_Data/SUMS_CompleteCatalog/Logs/"

# Test all step 1 files were ran 
if test1_step1files_completed:
    print("Testing all step 1 files were completed")
    for galaxy in ['lmc', 'smc']:
        path = f"H:/Data/SUMS_Tractor_Data/{galaxy}/"
        folders = glob.glob(f'{path}*X/*/')
        n = len(path) 
        save_names = []
        remaining = []
        for f in folders: 
            # File management
            obsid = f[n+15:n+20]
            uvfilter = f[n+23:n+26]
            segment = f[-4]
            extension = f[-2]
            step1_file_save = save_dir + f'Step1/{galaxy}/{obsid}_{uvfilter}_{segment}_{extension}_step1.csv'
            if os.path.exists(step1_file_save):
                continue
            else:
                remaining.append(f)
        if len(remaining) == 0:
            print(f'All files found for {galaxy}')
        else:
            print(f'{galaxy} has {len(remaining)/len(folders)*100:.2f}% remaining files to process')

if test2_step3_uniquecoords_areunique:
    print("Testing unique coordinates are unique")
    for galaxy in ['lmc', 'smc']:
        print(f"Testing {galaxy}")
        # Get files from step 2 
        files = glob.glob(step2_dir + f'{galaxy}_*.csv')
        for file_ in files:
            df = pd.read_csv(file_)
            coord = [(r,d) for r,d in zip(df.Ra,df.Dec)]
            unique_coord = np.array(list(dict.fromkeys(coord)))
            # Check that unique coords are unique 
            for co in unique_coord: 
                # Speed things up by getting a subset of the unique coordinates 
                threshhold = 0.001
                ra, dec = co
                ra_all, dec_all = unique_coord[:,0], unique_coord[:,1]
                sub_unique_coord = unique_coord[(np.abs(ra_all - ra) < threshhold) & (np.abs(dec_all - dec) < threshhold)]
                # Calculate distance between coordinates in arcseconds
                d = SkyCoord(co[0],co[1], unit='deg').separation(SkyCoord(sub_unique_coord[:,0],sub_unique_coord[:,1], unit='deg')).arcsecond
                # What unique coordinates are within X arcseconds?
                threshhold = 0.01
                min_d = d[d < threshhold]
                # Length above 1 to excluding self
                if len(min_d) > 1:
                    print(f'Multiple unique coordinates within {threshhold} arcseconds: {len(min_d)}')
                    print(sub_unique_coord[np.where(d < threshhold)])
                    print(min_d)
                    break
            print(f'Unique coordinates are unique in: {file_}')

if test3_step3_uniquecoords_assignment:
    print("Testing unique coordinates are assigned correctly")
    for galaxy in ['lmc', 'smc']:
        print(f"Testing {galaxy}")
        # Get files from step 2 
        files = glob.glob(step2_dir + f'{galaxy}_*.csv')
        for file_ in files:
            df = pd.read_csv(file_)
            df['group'] = 0
            counter = 1 
            # Will any coordinate be assigned different unique coordinates? 
            # If previously assigned, then the sum of the group numbers wouldn't be zero. 
            for coord in unique_coord:
                rows = df.loc[np.isclose(df.Ra,coord[0],atol=1e-5,rtol=1e-8) & np.isclose(df.Dec,coord[1],atol=1e-5,rtol=1e-8)]
                if np.sum(rows.group) > 0:
                    print("Multiple assignments: ", coord)
                    break
                df.loc[np.isclose(df.Ra,coord[0],atol=1e-5,rtol=1e-8) & np.isclose(df.Dec,coord[1],atol=1e-5,rtol=1e-8),'group'] = 1
                counter += 1
                
            # Will anything not be assigned in any coordinate? 
            # If after running there are still rows with group number 0,
            # then a unique coordinate was not assigned.
            unassigned = df[df.group == 0]
            if len(unassigned) > 0:
                print(f"{len(unassigned)} were not assigned a unique coordinate")
                break
            print(f'Unique coordinates are assigned correctly in {file_}')