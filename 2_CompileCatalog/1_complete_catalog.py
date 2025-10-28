import numpy as np 
import pandas as pd 
import time
import glob
import os 
from scipy.integrate import simpson
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
data_dir = os.getenv("DATADIR")
complete_catalog_dir = data_dir + "0_SUMS_Catalogs/CompleteCatalog/"

pd.options.mode.chained_assignment = None
# This code puts together a final catalog after Tractor and Heasarc have been run. 
# It is broken down into code blocks for ease of use.

############
## Inputs ##
############
galaxy = 'smc'

##############
## Switches ##
##############
# Drop NaNs, sources with high errors
# Calculate how crowded the source is
# - number of other sources within 5", 2.5", 1"
# - distance to the closest object
# - fraction of flux that was assigned to this source by tractor.
step_1 = False
# Combines CSVs with all the individual sources from each of the fields. 
# Adds a column that says which image a given measurement came from. 
step_2 = False
# Average together different observations of the same source to have a final measurement. 
# - Drop certain measurements:
#   -  Anything near an edge with an 'SSS' senstivity flag 
#   -  If it has a "residual fraction" greater than 0.3
# - Calculate the weighted average (weighted by residuals) 
# - Outputs: 
#   - the statistical error from tractor, propogated through the average 
#   - the standard deviation of the mean. 
#   - averages for the residuals
#   - fraction of the flux that was assigned to this source in each image.
# - Write out a new file, one line per source, with the new averaged photometry. 
step_3 = False
# Create the full catalog
# Link up the UVW1, UVM2, and UVW2 photometry for a given source. 
step_4 = True

###############
## Functions ##
###############

def flux_fraction_dis(distance,fwhm=2.5): # MRD
    #distance can be a list. In arcsec.

    # Generate the 2D gaussian (should be normalized):
    x1 = np.linspace(-5, 5, 41)
    y1 = np.linspace(-10, 10, 81)
    x, y = np.meshgrid(x1, y1)
    d_uvot = np.sqrt(x * x + y * y)
    index_out = np.where(d_uvot > 5.)
    
    frac = np.zeros(len(distance))
    for j,dist in enumerate(distance):

        d_trial = np.sqrt(x*x + (y - dist) * (y - dist))
        
        #Initialize the gaussian for the new source
        sigma = fwhm/2.355
        mu = 0
        g = (1./(2.*np.pi*sigma**2)) * np.exp(-((d_trial - mu) ** 2 / (2.0 * sigma ** 2)))
        
        #Zero out portions that aren't w/in 5" of the uvot position.
        g[index_out] = 0
    
        #Integrate the remaining flux.
        frac_x = [simpson(y=g[i],x=x1) for i in range(len(y1))]
        frac[j] = simpson(y=frac_x,x=y1)

    return frac
def weighted_avg_and_std(values, weights): # MRD
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return average, np.sqrt(variance)
def distance(x0,x1,y0,y1):
    return np.sqrt((x0-x1)**2+(y0-y1)**2)
def merge_with_tolerance(df1,df2,tolerance):
    temp_df1 = df1 / merge_with_tolerance
    
    temp_df1['ra'] = temp_df1['ra'].astype('int')
    temp_df1['dec'] = temp_df1['dec'].astype('int')


    temp_df2 = df2 / merge_with_tolerance
    temp_df2['ra'] = temp_df2['ra'].astype('int')
    temp_df2['dec'] = temp_df2['dec'].astype('int')

    merged = temp_df1.merge(temp_df2, how='outer', on=['ra','dec'])

    merged = merged * tolerance

    return merged
def mad(values,ra,dec,uvfilter,log_dir):
    """
    Returns the modified Z-score of all points in the array. This is defined as:
    Zm = 0.6745 * (xi - median(xi))/MAD
    where the MAD is the median average deviation defined as:
    MAD = median(abs(xi - median(xi)))
    """
    median = np.median(values)
    distances = values - median
    mad = np.median(abs(distances))

    # In rare cases it's possible that the values equal the median if more than half of the values are the same. 
    # This checks for that and if true returns values of 1 for zm rather than inf or nan.
    # I want to keep a text file of when this occurs. Shouldn't be often. 
    if mad != 0.0: 
        zm = 0.6745 * distances / mad 
        return np.abs(zm)
    else:
        filename = log_dir + 'no_mad_values_obtained.txt'
        line = f'ra: {ra} dec: {dec} uvfilter: {uvfilter}\n'
        if not os.path.exists(filename):
            file_object = open(filename, 'w+')
            file_object.write(line)
            file_object.close()
        else: 
            with open(filename, "a") as file_object:
                file_object.write(line)
        return np.ones(len(values))
# Only search close sources within coord_width degrees
def get_nearby_rows(row, df, width):
    rows = df[(df['ra']<row['ra']+width) & (df['ra']>row['ra']-width) & (df['dec']<row['dec']+width) & (df['dec']>row['dec']-width)]
    return rows
# If the closest source is within max_dist arcsec, return it
def get_closest_row(row, df):
    # Check all sources that are within a tenth of an arcsecond
    # Based on fiddling with topcat self cross match
    max_dist = 0.1 # arcsec
    coord_width = 0.0025 # degrees

    rows = get_nearby_rows(row, df, width=coord_width)
    # Remove self
    rows = rows[rows['temp_id'] != row['temp_id']]
    if len(rows) > 0:
        coords1 = SkyCoord(ra=row['ra'], dec=row['dec'], unit='deg')
        coords2 = SkyCoord(ra=rows['ra'].values, dec=rows['dec'].values, unit='deg')
        sep = coords1.separation(coords2).arcsec
        min_idx = np.argmin(sep)
        min_dist = np.min(sep)
        if min_dist <= max_dist:        
            return rows.iloc[min_idx], min_dist
        return None, None
    return None, None

############
## Step 1 ##
############

total_start = time.time()


if step_1:        
    path = f"H:/Data/SUMS_Tractor_Data/{galaxy}/"
    cat_path = f"H:/Data/SUMS_Tractor_Data/MatchedByPixel/{galaxy}/"
    folders = glob.glob(f'{path}*X/*/')
    n = len(path)

    counter = 0 
    for f in folders: 
        # File management
        obsid = f[n+15:n+20]
        uvfilter = f[n+23:n+26]
        segment = f[-4]
        extension = f[-2]
        step1_file_save = complete_catalog_dir + f'Step1/{galaxy}/{obsid}_{uvfilter}_{segment}_{extension}_step1.csv'

        # Code Start 
        if os.path.exists(step1_file_save):
            print('Skipping: ',step1_file_save)
            continue
        else:
            print(f"Running Step 1 on Obsid: {obsid} UVfilter:{uvfilter}")
            # Original image output from tractor 
            im_file = f + f'/{obsid}_{uvfilter}_{segment}_{extension}_img.fits'
            # Matched by pixel csv
            cat_file = cat_path + f'/{obsid}_{uvfilter}_{segment}_{extension}.csv'
            csv_data = pd.read_csv(cat_file)
            hdr = fits.open(im_file)[0]

            start = time.time()
            # Drop Nans
            csv_data = csv_data[~np.isnan(csv_data.MAG_ERR)].reset_index(drop=True)
            # Calculate how much tractor moved a source and add this as a row. 
            d = distance(csv_data.PIX_X,csv_data.INIT_PIX_X,csv_data.PIX_Y,csv_data.INIT_PIX_Y)
            csv_data['d_moved'] = d

            # Loop through sources
            for index, row in csv_data.iterrows():
                # If mag err > 3 sigma continue, don't calculate stats. 
                if csv_data.loc[index,'MAG_ERR'] > 0.36:
                    continue
                width = 20
                close_sources = csv_data[(row.INIT_PIX_X < csv_data.INIT_PIX_X + width) &
                                         (row.INIT_PIX_X > csv_data.INIT_PIX_X - width) &
                                         (row.INIT_PIX_Y < csv_data.INIT_PIX_Y + width) &
                                         (row.INIT_PIX_Y > csv_data.INIT_PIX_Y - width) ]

                # Drop the source itself from this subset 
                close_sources = close_sources.drop(index)

                d = distance(row.PIX_X,close_sources.PIX_X,row.PIX_Y,close_sources.PIX_Y)
               
                # Get all sources within 10" of the source, excluding the source itself.
                close_mag = close_sources.loc[(d< 10) & (d > 0),'MAG']

                if len(close_mag) > 0:
                    # Calculate flux from nearby sources within 10"
                    frac = flux_fraction_dis(d[(d< 10) & (d > 0)])
                    flux_from_others = 10.**(-0.4*close_mag) * frac
                    # Flux from the target
                    frac_self = flux_fraction_dis([0.0])
                    mag_self = csv_data.loc[index,'MAG']
                    flux_from_self =  10.**(-0.4*mag_self) * frac_self
                    # Calculate fraction of flux attributed to source.
                    flux_frac = flux_from_self/(flux_from_self + np.sum(flux_from_others))
                else:
                    # If there are no other sources within 20" then the fraction of flux is 1.0
                    flux_frac = 1.0 
                    
                # Store the number of sources around each source by distance
                csv_data.loc[index,'num1'] = len(close_sources[d < 1])
                csv_data.loc[index,'num2p5'] = len(close_sources[d < 2.5])
                csv_data.loc[index,'num5'] = len(close_sources[d < 5])

                # Store the closest source
                if len(d) > 0:
                    csv_data.loc[index,'closest'] = np.min(d)
                else:
                    csv_data.loc[index,'closest'] = np.nan

                # Store the fraction of flux attributed to the source. 
                csv_data.loc[index,'flux_frac'] = flux_frac

            # Drop Errors > 3 sigma
            csv_data = csv_data[csv_data.MAG_ERR < 0.36]
            
            # Drop if the source moved too much in an image (Like 1")
            csv_data = csv_data[csv_data['d_moved'] < 1]
            
            # Reset Index
            csv_data = csv_data.reset_index(drop=True)

            timetaken = (time.time() - start) / 60
            print(f'Time Taken: {timetaken} m')  
            csv_data.to_csv(step1_file_save,index=False)

            if counter == 0: 
                print('Estimated time to completion: ',timetaken * len(folders) / 60, ' hours')
            counter += 1

############
## Step 2 ##
############

if step_2:   
    path = f"H:/Data/SUMS_Tractor_Data/{galaxy}/"
    step1_dir = complete_catalog_dir  + "Step1/"
    step2_dir = complete_catalog_dir  + "Step2/"
    n = len(path)

    # Loop over each filter
    for uvfilter in ['um2','uw2','uw1']:
        step2_file_save = step2_dir + f'{galaxy}_{uvfilter}_step2.csv'
        if os.path.exists(step2_file_save):
            print('Skipping: ',step2_file_save)
            continue

        # Get all the folders with data from that filter
        folders = glob.glob(f'{path}*X/*{uvfilter}*/')

        start = time.time()
        csvs = []
 
        for f in folders: 
            # File management
            obsid = f[n+15:n+20]
            uvfilter = f[n+23:n+26]
            segment = f[-4]
            extension = f[-2]
            im_file = f + f'/{obsid}_{uvfilter}_{segment}_{extension}_img.fits'
            step1_cat_file = step1_dir + f'{galaxy}/{obsid}_{uvfilter}_{segment}_{extension}_step1.csv'

            # Open Step 1 file 
            csv_data = pd.read_csv(step1_cat_file)
            hdr = fits.open(im_file)[0]
            # Get exposure time
            exp_time = hdr.header["EXPOSURE"]
            
            print(f"Running Step 2 on Obsid: {obsid} UVfilter:{uvfilter}")

            # Add image name and exposure as a column 
            csv_data['filename'] = np.repeat(f'{obsid}_{segment}_{extension}',csv_data.shape[0])
            csv_data['exposure'] = np.repeat(exp_time,csv_data.shape[0])
        
            csvs.append(csv_data)

        combined_csv_data = pd.concat(csvs)
        combined_csv_data.to_csv(step2_file_save,index=False)
        print(f'Time taken for {uvfilter}: {time.time()-start} s')

############
## Step 3 ##
############

if step_3: 
    step2_dir = complete_catalog_dir  + "Step2/"
    step3_dir = complete_catalog_dir + "Step3/"
    log_dir = data_dir + "0_SUMS_Catalogs/Logs/"

    # Get files from step 2 
    s2_um2 = pd.read_csv(step2_dir + f'{galaxy}_um2_step2.csv')
    s2_uw1 = pd.read_csv(step2_dir + f'{galaxy}_uw1_step2.csv')
    s2_uw2 = pd.read_csv(step2_dir + f'{galaxy}_uw2_step2.csv')

    # From 0_CleanData get list of files to skip
    bad_files = pd.read_csv('../0_CleanData/files_to_skip.csv')

    bad_files = bad_files[bad_files.galaxy == galaxy]

    # Drop images with tracking issues 
    s2_um2 = s2_um2[~s2_um2.filename.isin(bad_files.loc[bad_files['filter'] == 'uvm2', 'files'])].reset_index(drop=True)
    s2_uw1 = s2_uw1[~s2_uw1.filename.isin(bad_files.loc[bad_files['filter'] == 'uvw1', 'files'])].reset_index(drop=True)
    s2_uw2 = s2_uw2[~s2_uw2.filename.isin(bad_files.loc[bad_files['filter'] == 'uvw2', 'files'])].reset_index(drop=True)
    
    files = [s2_um2,s2_uw1,s2_uw2]
    uvfilters = ['uvm2','uvw1','uvw2']

    # Loop over each file and do some reductions 
    for df,uvfilter in zip(files,uvfilters): 
        print(f'Running step 3 for uvfilter: {uvfilter}')
        savename = step3_dir+ f'{galaxy}_{uvfilter}_step3.csv'
        # Drop bad rows / keep good rows, save how much we're dropping  
        s2_shape = df.shape
        # Save here so we can see what the RESID_Frac cut will drop
        df.to_csv(log_dir+f'{galaxy}_{uvfilter}_resid_frac_TOTAL.csv',index=False)
        # Otherwise drop
        df = df[np.abs(df.RESID_FRAC) < 0.3]
        rf_shape = df.shape
        df =df[df.SSS == 1.0]
        sss_shape = df.shape
        df =df[df.EDGE == 1.0]
        edge_shape = df.shape

        # Write out how much we're dropping 
        phrases = [f'Reductions for {uvfilter} \nShape after step 2: {s2_shape}',
                   f'Shape after dropping resid_frac > 0.3: {rf_shape}',
                   f'Shape after dropping sss != 1: {sss_shape}',
                   f'Shape after dropping edge!= 1: {edge_shape}']

        with open(log_dir+f'reductions_log_for_{galaxy}_{uvfilter}.txt','w') as f:
            for phrase in phrases:
                f.write(phrase)
                f.write('\n')

        # Get all unique coordinates
        coord = [(r,d) for r,d in zip(df.Ra,df.Dec) ] 
        coord_unique = np.array(list(dict.fromkeys(coord)))
        # Counter
        counter = 0
        # For each unique coordinate get a group of that coordinate
        rows = []
        for r,d in coord_unique:
            # Get Group
            g = df[np.isclose(df.Ra,r,atol=1e-5,rtol=1e-8) & np.isclose(df.Dec,d,atol=1e-5,rtol=1e-8)]
            length = len(g)
            # Start a dictionary:
            # Combine n neighbors 
            # MCPS should be all the same for each member of the group so take the first. 
            line = {'ra':r,
                    'dec':d,
                    'U':g.Umag.values[0],
                    'e_U':g.e_Umag.values[0],
                    'B':g.Bmag.values[0],
                    'e_B':g.e_Bmag.values[0],
                    'V':g.Vmag.values[0],
                    'e_V':g.e_Vmag.values[0],
                    'I':g.Imag.values[0],
                    'e_I':g.e_Imag.values[0],
                    'J':g.Jmag.values[0],
                    'e_J':g.e_Jmag.values[0],
                    'H':g.Hmag.values[0],
                    'e_H':g.e_Hmag.values[0],
                    'Ks':g.Ksmag.values[0],
                    'e_Ks':g.e_Ksmag.values[0]}
            if length == 1: 
                    line['num_obs'] = 1.
                    line['num_outliers'] = 0.
                    line['mag'] = g.MAG.values[0]
                    line['mag_err'] = g.MAG_ERR.values[0]
                    line['flux_frac'] = g.flux_frac.values[0]
                    line['resid_frac'] = g.RESID_FRAC.values[0]

                    # Stuff that doesn't apply if there's only one observation. 
                    line['mag_std'] = np.nan
                    line['flux_frac_std'] = np.nan
                    line['resid_frac_std'] = np.nan
                    line['std_unweighted'] = np.nan

            # If there is more than one line in the group
            else:
                # Get the MAD Zm values
                zm = mad(g.MAG.values,r,d,uvfilter,log_dir)

                # Reject anything above 3.5 / Keep anything below 3.5
                n_rejected = len(g[zm > 3.5])
                g = g[zm < 3.5]

                # Calculate weights and magnitudes 
                errors = np.sqrt(g.RESID_FRAC**2 + g.MAG_ERR**2)
                weights = 1./(errors**2)
                mean_mag, std_mag = weighted_avg_and_std(g.MAG,weights)

                # Calculate weighted fluxfrac residfrac and magerror
                mean_flux_frac, std_flux_frac = weighted_avg_and_std(g.flux_frac,weights)
                mean_resid, std_resid = weighted_avg_and_std(g.RESID_FRAC,weights)

                # Calculate magnitude errors 
                mean_mag_error= np.sqrt(np.sum(weights**2*g.MAG_ERR**2))/np.sum(weights) 
                
                line['num_obs'] = length
                line['num_outliers'] = n_rejected
                line['mag'] = mean_mag
                line['mag_err'] = mean_mag_error
                line['mag_std'] = std_mag
                line['flux_frac'] = mean_flux_frac
                line['flux_frac_std'] = std_flux_frac
                line['resid_frac'] = mean_resid
                line['resid_frac_std'] = std_resid
                line['std_unweighted'] = np.std(g.MAG)

            # How many neighbors could be around?             
            line['num5'] = max(g.num5) 
            line['num2p5'] = max(g.num2p5)
            line['num1'] = max(g.num1)
            # What is the closest neighbors?
            line['closest_min'] = np.min(g.closest)
            # How far does tractor typically have to move this source? 
            line['dist_moved'] = np.mean(g.d_moved)
            line = pd.DataFrame([line])
            # For first line start a new csv, after that append each line to it. 
            rows.append(line)
            counter += 1
            if counter % 1000 == 0:
                print (f'{counter / len(coord_unique) * 100:.2f} % Complete')

        reduced_csv = pd.concat(rows)
        reduced_csv.to_csv(savename,index=False)
        print(f'Finished {galaxy} {uvfilter}.') 


############
## Step 4 ##
############
if step_4:
    step3_dir = complete_catalog_dir  + "Step3/"
    step4_dir = complete_catalog_dir  + "Step4/"

    print(f'Running step 4 on {galaxy}')

    #### FORMAT COLUMNS AND MERGE ####
    filter_specific_keys = ['num_obs', 'num_outliers','mag', 'mag_err', 'flux_frac', 'resid_frac', 'mag_std', 'flux_frac_std', 'resid_frac_std', 'std_unweighted', 
                            'num5', 'num2p5', 'num1', 'closest_mean', 'closest_min', 'closest_std', 'dist_moved']

    # Rename each filter specific key to have the filter in front of it 
    def rename_keys(df,filter_name):
        for key in filter_specific_keys:
            df = df.rename(columns={key:f'{filter_name}_{key}'})
        return df
    
    uv_filters = ['uvw2','uvm2','uvw1']
    dfs = [] 
    mcps_cols = ['ra','dec','U','e_U','B','e_B','V','e_V','I','e_I','J','e_J','H','e_H','Ks','e_Ks']
    mcps_drop_cols =  ['J','e_J','H','e_H','Ks','e_Ks']
    # Loop over each filter
    for uvfilter in uv_filters:
        # Read in the step 3 files, ignore unnamed column by setting index_col=0 and resetting the index 
        df = pd.read_csv(step3_dir+f'{galaxy}_{uvfilter}_step3.csv')
        # Rename the keys to have the uvfilter in front 
        df = rename_keys(df,uvfilter)
        # Turn any spaces to NaNs to make consistent
        df = df.mask(df=='')
        # Turny any mcps filters that are equal to 0 or -99 to NaN
        df[mcps_cols] = df[mcps_cols].replace({0:np.nan,-99:np.nan})
        # Save the file
        dfs.append(df)

    # Combine the dataframes together
    # How = outer, means dont delete anything https://stackoverflow.com/questions/53645882/pandas-merging-101
    step_4_df = pd.merge(dfs[0],dfs[1],on=mcps_cols,how='outer').merge(dfs[2],on=mcps_cols,how='outer')
    # Drop some mcps cols 
    step_4_df = step_4_df.drop(columns=mcps_drop_cols)
    step_4_df = step_4_df.reset_index(drop=True)


    #### DROP DUPLICATES ####
    step_4_df['temp_id'] = np.arange(len(step_4_df))

    new_df = step_4_df.copy()
    # Get a list of columns ignoring ra/dec because of float precision and the new key we just made
    columns = step_4_df.columns[2:-1]
    uv_columns = columns[8:]
    optical_columns = columns[:8]

    keep = []
    drop = []
    for ind, row in step_4_df.iterrows():
        # Progress print
        if ind % 100_000 == 0:
            print(f"Processed {ind/step_4_df.shape[0]*100:.2f}% of rows")

        # Get the closest row
        closest_row, dist = get_closest_row(row, step_4_df)

        dupe_id = int(closest_row['temp_id']) if closest_row is not None else None

        # No close row within threshold
        if dupe_id is None:
            continue
            
        # Did we already account for it?
        if (dupe_id in keep) or (dupe_id in drop):
            continue
        
        # If rows are exactly the same then we know theyre duplicates
        same_rows = step_4_df.loc[ind, columns].equals(step_4_df.loc[dupe_id, columns])
        if same_rows == False:
            print(f"Rows are completely identical")

        # If they are not exactly the same, they may have some missing UV values
        # First check if the optical values are the same, 
        # then replace nan values within finite value in other row

        # MCPS has its own duplicates, but we don't remove those here.
        else:
            # If the optical values are not the same but they're both finite values then don't consider it a duplicate
            for col in optical_columns:
                if (step_4_df.loc[ind, col] != step_4_df.loc[dupe_id, col]) and np.isfinite(step_4_df.loc[ind, col]) and np.isfinite(step_4_df.loc[dupe_id, col]):
                    continue

            # If the dupe row has values the other row doesnt, fill them in
            for col in uv_columns:
                if pd.isna(step_4_df.loc[ind, col]) and not pd.isna(step_4_df.loc[dupe_id, col]):
                    new_df.loc[ind, col] = step_4_df.loc[dupe_id, col]
            
            # Check theyre the same now, including if both are nans
            same_rows = new_df.loc[ind, columns].equals(new_df.loc[dupe_id, columns])
            
            # The above line doesn't always work, may be a data type issue, so double check each column
            if same_rows == False:
                # First check if the different rows are all nans in one of the rows
                for col in uv_columns:
                    if new_df.loc[ind,col] != new_df.loc[dupe_id,col]:
                        # Are both finite?
                        if np.isfinite(new_df.loc[ind,col]) and np.isfinite(new_df.loc[dupe_id,col]):
                            print(f"After filling, row {ind} is still not identical to row {dupe_id} for column {col}")

        # Keep the current row for comparison when testing
        keep.append(ind)
        # Drop the other row
        drop.append(dupe_id)

    print(f"Dropping {len(drop)} duplicate rows out of {step_4_df.shape[0]} total rows")
    new_df = new_df[~new_df['temp_id'].isin(drop)].reset_index(drop=True)

    # Drop temp_id column
    new_df = new_df.drop(columns=['temp_id'])

    # Save
    new_df.to_csv(step4_dir+f'{galaxy}_photometry.csv',index=False)

##########
# Finish # 
##########
print(f'Total Time Taken: {(time.time() - total_start)/60} minutes')
print("Step Complete") 
