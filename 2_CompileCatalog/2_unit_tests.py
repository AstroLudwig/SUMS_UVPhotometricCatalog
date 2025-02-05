import glob 
import os
import pandas as pd 
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits

step1_files_completed = False
step1_fileshavecorrectamountofdata = False
step2_fileshavecorrectamountofdata = False
step3_uniquecoords_areunique = False
step3_uniquecoords_assignment = True

# Files 
save_dir = "C:/Projects/0_Data/0_SUMS_Catalogs/CompleteCatalog/"
step2_dir = "C:/Projects/0_Data/0_SUMS_Catalogs/CompleteCatalog/Step2/"
log_dir = "C:/Projects/0_Data/0_SUMS_Catalogs/CompleteCatalog/Logs/"

# Test all step 1 files were ran 
if step1_files_completed:
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
                # Open it and make sure theres data
                try:
                    df = pd.read_csv(step1_file_save)
                    continue
                except:
                    print(f'{step1_file_save} has no data')
                    break
            else:
                remaining.append(f)
        if len(remaining) == 0:
            print(f'All files found for {galaxy}')
        else:
            print(f'{galaxy} has {len(remaining)/len(folders)*100:.2f}% remaining files to process')

if step1_fileshavecorrectamountofdata:   
    print("Testing all step 1 files have correct amount of data")
    def distance(x0,x1,y0,y1):
        return np.sqrt((x0-x1)**2+(y0-y1)**2)
    for galaxy in ['lmc', 'smc']:  
        path = f"H:/Data/SUMS_Tractor_Data/{galaxy}/"
        cat_path = f"H:/Data/SUMS_Tractor_Data/MatchedByPixel/{galaxy}/"
        save_dir = "C:/Projects/0_Data/0_SUMS_Catalogs/CompleteCatalog/"
        folders = glob.glob(f'{path}*X/*/')
        n = len(path) 
        counter = 0 
        for f in folders: 
            # File management
            obsid = f[n+15:n+20]
            uvfilter = f[n+23:n+26]
            segment = f[-4]
            extension = f[-2]
            step1_file_save = save_dir + f'Step1/{galaxy}/{obsid}_{uvfilter}_{segment}_{extension}_step1.csv'

            # Original image output from tractor 
            im_file = f + f'/{obsid}_{uvfilter}_{segment}_{extension}_img.fits'
            # Matched by pixel csv
            cat_file = cat_path + f'/{obsid}_{uvfilter}_{segment}_{extension}.csv'
            csv_data = pd.read_csv(cat_file)
            hdr = fits.open(im_file)[0]

            # Drop Nans
            csv_data = csv_data[~np.isnan(csv_data.MAG_ERR)].reset_index(drop=True)
            # Calculate how much tractor moved a source and add this as a row. 
            d = distance(csv_data.PIX_X,csv_data.INIT_PIX_X,csv_data.PIX_Y,csv_data.INIT_PIX_Y)
            csv_data['d_moved'] = d

            # Drop Errors > 3 sigma
            csv_data = csv_data[csv_data.MAG_ERR < 0.36]
            
            # Drop if the source moved too much in an image (Like 1")
            csv_data = csv_data[csv_data['d_moved'] < 1]
            
            # Reset Index
            csv_data = csv_data.reset_index(drop=True)

            # Open the saved file from step 1 
            df = pd.read_csv(step1_file_save)

            # Compare rows 
            if df.shape[0] != csv_data.shape[0]:
                print('Shapes do not match.',step1_file_save)
                break

if step2_fileshavecorrectamountofdata:
    for galaxy in ['lmc', 'smc']:
        path = f"H:/Data/SUMS_Tractor_Data/{galaxy}/"
        step1_dir = "C:/Projects/0_Data/0_SUMS_Catalogs/CompleteCatalog/Step1/"
        step2_dir = "C:/Projects/0_Data/0_SUMS_Catalogs/CompleteCatalog/Step2/"
        n = len(path)

        # Loop over each filter
        for uvfilter in ['um2','uw2','uw1']:
            step2_file_save = step2_dir + f'{galaxy}_{uvfilter}_step2.csv'

            # Get all the folders with data from that filter
            folders = glob.glob(f'{path}*X/*{uvfilter}*/')

            # Open all files and get the number of rows
            n_rows = []
            for f in folders: 
                # File management
                obsid = f[n+15:n+20]
                uvfilter = f[n+23:n+26]
                segment = f[-4]
                extension = f[-2]
                im_file = f + f'/{obsid}_{uvfilter}_{segment}_{extension}_img.fits'
                step1_cat_file = step1_dir + f'{galaxy}/{obsid}_{uvfilter}_{segment}_{extension}_step1.csv'
                csv_data = pd.read_csv(step1_cat_file)
                n_rows.append(csv_data.shape[0])

            # Open the saved file from step 2
            df = pd.read_csv(step2_file_save)
            # Compare to make sure they match
            if df.shape[0] != sum(n_rows):
                print('Shapes do not match.',step2_file_save)
                break
        print(f'{galaxy} has the correct amount of data')

if step3_uniquecoords_areunique:
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

if step3_uniquecoords_assignment:
    print("Testing unique coordinates are assigned correctly")
    for galaxy in ['lmc', 'smc']:
        print(f"Testing {galaxy}")
        # Get files from step 2 
        files = glob.glob(step2_dir + f'{galaxy}_*.csv')
        for file_ in files:
            df = pd.read_csv(file_)
            coord = [(r,d) for r,d in zip(df.Ra,df.Dec)]
            unique_coord = np.array(list(dict.fromkeys(coord)))
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