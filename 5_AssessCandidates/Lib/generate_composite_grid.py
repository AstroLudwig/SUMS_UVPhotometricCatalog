import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob

# This file creates a composite grid of MIST isochrones with Ylva's stripped star models
# Models Accessed with: https://waps.cfa.harvard.edu/MIST/interp_tracks.html
# For metallicities: 
# LMC-like: log10(0.006/0.0142) = -0.37
# SMC-like: log10(0.002/0.0142) = -0.85
# Masses: 
# [1, 12.4, steps=0.2]

class CompositePhotometry:

    def __init__(self,galaxy,data_dir,save_csv=None,rotation='Rotation'):

        # Galaxy specific info 
        self.galaxy = galaxy.upper()
        self.synth_cols = ['Minit_strip', 'M_strip', 'M_MS', 'frac_MS']
        self.phot_filters = ['UVW2','UVM2','UVW1','U','B','V','I']
      
        if self.galaxy == 'LMC':
            self.distance = 50e3
            # Synthetic Photometry 
            self.synth_photometry_file = data_dir+'1_Models/StrippedStars/photometry_CMFGEN_composites.csv'
            self.synth_photometry = pd.read_csv(self.synth_photometry_file)
            self.synth_photometry.rename(columns={col.lower():col.upper() for col in self.phot_filters},inplace=True)
            # Stripped Star Evolutionary Model - For getting stripped star mass from initial mass
            self.stripped_evolution_file = data_dir+'1_Models/Gotberg18/Gotberg2018Z_0p006_LMC_properties.csv'
            self.synth_photometry = self.Handle_Synth_Photometry(to_apparent=False)
            # For the isolated stripped stars
            self.synth_stripped = self.synth_photometry[self.synth_photometry['M_MS'] == 0].copy()
            
        elif self.galaxy == 'SMC':
            self.distance = 60.6e3
            self.synth_photometry_file = data_dir+'1_Models/StrippedStars/stripped_stars_Z0.002_ABmag.txt'
            self.synth_photometry = pd.read_csv(self.synth_photometry_file, comment='#',delimiter='\\s+')
            # Rename Minit to 'Minit_strip'
            self.synth_photometry.rename(columns={'Minit':'Minit_strip'},inplace=True)
            self.synth_photometry['frac_MS'] = 0.0
            self.synth_photometry['M_MS'] = 0.0
            self.stripped_evolution_file = data_dir+'1_Models/Gotberg18/Gotberg2018Z_0p002_SMC_properties.csv'
            self.synth_photometry = self.Handle_Synth_Photometry()
            # For the isolated stripped stars
            self.synth_stripped = self.synth_photometry.copy()

        # Mist zero points for converting from Vega to AB
        self.vega_ab_table = pd.read_csv('https://waps.cfa.harvard.edu/MIST/BC_tables/zeropoints.txt',delimiter='\\s+')
        
        # MIST file names for Swift Photometry
        self.rotation = rotation
        self.mist_files = glob.glob(f'{data_dir}/1_Models/MIST/{self.galaxy}/MIST/{self.rotation}/*.track.*')     

        # Main Sequence Photometry
        self.ms_df = self.Generate_MS_Table(data_dir)

        # Composite Table
        self.composite_df = self.Generate_Composite_Table()

        # Save the composite table
        if save_csv:
            self.composite_df.to_csv(save_csv,index=False)

    def convert_to_flux(self,mag):
        # Convert from AB mag to flux (For Getting Fractional FLux Values)
        return 10**(-1*(mag+48.6)/2.5)
    
    def get_fractional_flux(self,stripped_mag,ms_mag):
        combined_mag = self.Combine_Mags(stripped_mag,ms_mag)
        stripped_flux = self.convert_to_flux(stripped_mag)
        combined_flux = self.convert_to_flux(combined_mag)
        return stripped_flux/combined_flux
    
    def Generate_Composite_Table(self):
        # Create a grid of all possible MS and Stripped Star masses       
        M_MS = self.ms_df.M_MS.values
        M_SHS = self.synth_stripped.M_strip.values
        M_MS, M_SHS = np.meshgrid(M_MS, M_SHS)
        M_MS = M_MS.flatten()
        M_SHS = M_SHS.flatten()

        rows = []
        for m_ms, m_shs in zip(M_MS,M_SHS):
            # Get the MS and SHS row
            synth_row = self.synth_stripped.loc[self.synth_stripped.M_strip == m_shs]
            mist_row = self.ms_df.loc[self.ms_df.M_MS == m_ms]
            # Store relevant info and combined magnitudes
            row = {'Minit_strip':synth_row.Minit_strip.values[0],
                'M_strip':m_shs,
                'M_MS':np.round(m_ms,2),
                'frac_MS':0.2,
                'UVW2': self.Combine_Mags(mist_row.UVW2.values[0],synth_row.UVW2.values[0]),
                'UVM2': self.Combine_Mags(mist_row.UVM2.values[0],synth_row.UVM2.values[0]),
                'UVW1': self.Combine_Mags(mist_row.UVW1.values[0],synth_row.UVW1.values[0]),
                'U': self.Combine_Mags(mist_row.U.values[0],synth_row.U.values[0]),
                'B': self.Combine_Mags(mist_row.B.values[0],synth_row.B.values[0]),
                'V': self.Combine_Mags(mist_row.V.values[0],synth_row.V.values[0]),
                'I': self.Combine_Mags(mist_row.I.values[0],synth_row.I.values[0]),
                'UVW2_flux_frac':self.get_fractional_flux(synth_row.UVW2.values[0],mist_row.UVW2.values[0]),
                'UVM2_flux_frac':self.get_fractional_flux(synth_row.UVM2.values[0],mist_row.UVM2.values[0]),
                'UVW1_flux_frac':self.get_fractional_flux(synth_row.UVW1.values[0],mist_row.UVW1.values[0]),
                'U_flux_frac':self.get_fractional_flux(synth_row.U.values[0],mist_row.U.values[0]),
                'B_flux_frac':self.get_fractional_flux(synth_row.B.values[0],mist_row.B.values[0]),
                'V_flux_frac':self.get_fractional_flux(synth_row.V.values[0],mist_row.V.values[0]),
                'I_flux_frac':self.get_fractional_flux(synth_row.I.values[0],mist_row.I.values[0])}
            rows.append(row)

        # Create dataframe
        composite_df = pd.DataFrame(rows)

        # We also want cases have frac_MS = 0.0 (stripped star dominated) 
        # and cases where frac_MS = 1.0 (no stripped star)
        frac_columns = ['UVW2_flux_frac','UVM2_flux_frac','UVW1_flux_frac','U_flux_frac','B_flux_frac','V_flux_frac','I_flux_frac']
        for col in frac_columns:
            self.synth_stripped[col] = 1.0
            self.ms_df[col] = 0.0
        self.ms_df['Minit_strip'] = 0.0
        self.ms_df['M_strip'] = 0.0
        self.ms_df['frac_MS'] = 0.2

        # Combine synth stripped and composite grid so we also have frac_MS = 0.0
        composite_df = pd.concat([composite_df,self.synth_stripped])

        # Combine mist models and composite grid so we also have frac_MS = 1.0
        composite_df = pd.concat([composite_df,self.ms_df])

        # Sort by frac_MS, then Minit_strip, then M_MS
        composite_df = composite_df.sort_values(by=['frac_MS','Minit_strip','M_MS'])

        return composite_df

    def Handle_Synth_Photometry(self,to_apparent=True):        
        # Convert synthetic photometry for stripped stars to apparent magnitudes
        if to_apparent:
            for col in self.phot_filters:
                    self.synth_photometry[col] = self.Absolute_to_Apparent(self.synth_photometry[col])

        # Add the mass of the stripped star to the synthetic photometry
        self.stripped_evolution = pd.read_csv(self.stripped_evolution_file)
        for ind, row in self.stripped_evolution.iterrows():
            self.synth_photometry.loc[self.synth_photometry.Minit_strip.astype(str) == str(row["Minit"]), "M_strip"] = row["Mstrip"]
        
        # Rearrange the columns 
        return self.synth_photometry[self.synth_cols + self.phot_filters]

    def Generate_MS_Table(self,data_dir):
        rows = []
        for mist_file in self.mist_files:
            # Get the associated UBVI filename
            track = mist_file.split('/')[-1].split('.')[0]
            swift_file = glob.glob(f'{data_dir}1_Models/MIST/{self.galaxy}/Swift/{self.rotation}/{track}*.cmd')[0]
            ubvi_file = glob.glob(f'{data_dir}1_Models/MIST/{self.galaxy}/UBVI/{self.rotation}/{track}*.cmd')[0]
   
            # Get MIST data
            mist = self.Open_File(mist_file,11,mist=True)
            swift = self.Open_File(swift_file,14,mist=True)
            ubvi = self.Open_File(ubvi_file,14,mist=True)
            # Get the row closest to the frac_MS age
            age, frac_ms = self.Get_Age(mist)
            mist = self.Get_Row_By_Age(mist,age)
            swift = self.Get_Row_By_Age(swift,age)
            ubvi = self.Get_Row_By_Age(ubvi,age)

            # Check that it's the same age
            if mist.star_age != swift.star_age or mist.star_age != ubvi.star_age:
                print('Ages are not the same')

            # Get the row
            row = {'frac_MS':frac_ms,
                'M_MS':np.round(mist.star_mass,2),
                # Swift is already in AB
                'UVW2':swift.Swift_UVW2,
                'UVM2':swift.Swift_UVM2,
                'UVW1':swift.Swift_UVW1,
                # Convert from Vega to AB for UBVI
                'U':self.Vega2AB(ubvi,'Bessell_U'),
                'B':self.Vega2AB(ubvi,'Bessell_B'),
                'V':self.Vega2AB(ubvi,'Bessell_V'),
                'I':self.Vega2AB(ubvi,'Bessell_I')}

            rows.append(row)

        # Combine main sequence photometry into a dataframe
        ms_df = pd.DataFrame(rows)

        # Put all magnitudes into apparent magnitudes
        for col in ms_df.columns[2:]:
            ms_df[col] = self.Absolute_to_Apparent(ms_df[col])
        
        # Sort by mass 
        ms_df = ms_df.sort_values(by=['M_MS']).reset_index(drop=True)

        return ms_df


    # Read in files
    def Open_File(self,file_name,header_number,mist=False):
        # Get the column names
        columns = open(file_name, 'r').readlines()[header_number].split()[1:]
        # Check that it's correct
        if mist:
            if columns[0] != 'star_age':
                print('Header number is wrong')
        # Read in the file
        df = pd.read_csv(file_name, comment='#',names = columns,delimiter='\\s+')
        return df
    
    # Get the age of the MS star at some percent of it's total lifetime
    def Get_Age(self,mist,frac_MS=0.2):
        # Phase = 0 corresponds to the main sequence
        min_age = mist.loc[mist['phase'] == 0.0,'star_age'].min()
        max_age = mist.loc[mist['phase'] == 0.0,'star_age'].max()
        age = frac_MS * (max_age - min_age) + min_age
        frac_ms = np.round(age/max_age,2)
        return age, frac_ms
    
    def Get_Row_By_Age(self,df,age):
        # Find where frac_MS = 0.2
        age_diff = np.abs(df['star_age'] - age)
        row =  df.iloc[np.argmin(age_diff)]
        how_close = np.abs(row.star_age - age)/age * 100
        if how_close > 2:
            print(f'The closest age is {how_close}% away from age')
        return row
    
    def Vega2AB(self,df,band):
        # Swift is already in AB but MIST needs to be converted from Vega
        vega_ab = self.vega_ab_table.loc[self.vega_ab_table['filter'] == band,'mag(Vega/AB)']
        converted_mag = df[band] + vega_ab.values[0]
        return converted_mag
    
    def Absolute_to_Apparent(self,AbsoluteMag):
        # Convert Absolute to Apparent Magnitude
        #print(AbsoluteMag.values[0],self.distance)
        return AbsoluteMag + 5 * (np.log10(self.distance/10))
    
    
    def Combine_Mags(self,mag1,mag2):
        # Combine two magnitudes
        # These should be in apparent magnitudes
        return -2.5*np.log10(10**(-0.4*mag1) + 10**(-0.4*mag2))