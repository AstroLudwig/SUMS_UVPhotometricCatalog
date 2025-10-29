import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import numpy as np 
import os 
import matplotlib
from astropy.io import ascii
from matplotlib import rc
from scipy import stats
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy import interpolate
import os

# Various functions to speed up plotting

data_dir = os.getenv("DATADIR")
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.family'] = "serif"

class CMDFormatter:
    
    def __init__(self,galaxy,cmap_vals=[0.55,0.95],create_kde=False,create_zams=False,create_models=True):
        self.galaxy = galaxy
        
        # Get galaxy specific variables
        if self.galaxy == 'smc' or self.galaxy == 'SMC':
            self.get_smc(data_dir)
            
        elif self.galaxy == 'lmc' or self.galaxy == 'LMC':
            self.get_lmc(data_dir)
            

        # Create or Import ZAMS csv
        if create_zams:
            self.zams = self.create_zams_csv(data_dir)
        else:
            self.zams = pd.read_csv(data_dir + f'1_Models/ZAMS/{self.galaxy}_zams_apparent.csv')

        # Create or Import KDE csvs
        if create_kde:
            self.create_kde_files(data_dir)
        else:
            self.kdes = {'uvm2_v': pd.read_csv(data_dir+f'0_SUMS_Catalogs/KDEs/{galaxy}_uvm2_v_kde.csv'),
                         'uvw2_u': pd.read_csv(data_dir+f'0_SUMS_Catalogs/KDEs/{galaxy}_uvw2_u_kde.csv'),
                         'uvw1_b': pd.read_csv(data_dir+f'0_SUMS_Catalogs/KDEs/{galaxy}_uvw1_b_kde.csv'),
                         'u_v': pd.read_csv(data_dir+f'0_SUMS_Catalogs/KDEs/{galaxy}_u_v_kde.csv')}

        # Composite Grid 
        self.composites = pd.read_csv(data_dir+f'0_SUMS_Catalogs/CompositeGrid/{galaxy.lower()}_composite_photometry.csv')
        self.format_composites()

        # Create KDE colormap
        self.kde_cmap = self.truncate_colormap(plt.get_cmap('Greys'),cmap_vals[0], cmap_vals[1])

        # SUMS Candidates
        self.candidates = pd.read_csv(data_dir+'0_SUMS_Catalogs/CandidateCatalog/3_stripped_star_candidates.csv')
        self.format_candidates()

        # Mist Models
        self.mist_dir =  data_dir+f'1_Models/MIST/{galaxy.upper()}/'


        # Which masses for evolutionary models to show 
        self.mist_fnames = ['0030000M','0036000M','0042000M','0048000M', '0054000M','0060000M','0066000M','0072000M','0078000M','0084000M',
                            '0090000M','0096000M','0102000M','0108000M','0114000M','0120000M']

    #############################
    # Galaxy specific variables #
    #############################
    def get_smc(self,data_dir):
        self.zams_file = 'ZAMS_Z0.002_ABmag.txt'
        self.distance = 60.6e3
        self.av = 0.22
        self.extinction_coeff = {
            'uvw2' : 3.522767303652845,
            'uvm2' : 3.062652780996107,
            'uvw1' : 2.7436011426496876,
            'u'  : 1.746851127566682,
            'b'  : 1.4073996800012445,
            'v'  : 1.0353852912271932,
            'i'  : 0.5911705440475679,
        }
        self.wr = pd.read_csv(data_dir+'7_SourceCatalogs/WR/SMC_WR_DG23.csv')
    def get_lmc(self,data_dir):
        self.zams_file = 'ZAMS_Z0.006_ABmag.txt'
        self.distance = 50e3
        self.av = 0.38
        self.extinction_coeff = {
            'uvw2' : 2.644541680555541,
            'uvm2' : 2.7233599880955124,
            'uvw1' : 2.3536449902431045,
            'u'  : 1.5634790597438197,
            'b'  : 1.3170082045312625,
            'v'  : 1.0333402844940784,
            'i'  : 0.7366865234305091,
        }
        self.wr = pd.read_csv(data_dir+'7_SourceCatalogs/WR/LMC_WR_2024.csv')
    ############################
    # Absolute to Apparant Mag #
    ############################
    def Absolute_to_Apparent(self,AbsoluteMag,distance=None):
        if distance == None:
            distance = self.distance
        return AbsoluteMag + 5 * (np.log10(distance/10))

    #########################
    # Extinction Correction #
    #########################
    def deredden(self,mag,band):
        return mag - self.extinction_coeff[band] * self.av
    
    ################################################
    # Create csv of ZAMS in apparent AB magnitudes #
    ################################################
    def create_zams_csv(self,data_dir):
        # Import Zams Model 
        zams = pd.read_csv(data_dir + f'1_Models/ZAMS/{self.zams_file}',sep='\\s+',comment='#')
        # Get the columns we need 
        cols = ['Minit','UVW2_spec','UVM2_spec','UVW1_spec','U_spec','B_spec','V_spec','I_spec']
        zams = zams[cols]
        # Convert to apparent magnitudes
        for col in cols[1:]:
            zams[col] = self.Absolute_to_Apparent(zams[col])
        # Rename columns
        new_cols = ['m','uvw2','uvm2','uvw1','u','b','v','i']
        zams.columns = new_cols
        # Save as dataframe
        zams.to_csv(data_dir + f'1_Models/ZAMS/{self.galaxy}_zams_apparent.csv')
        return zams

    #################
    # KDE Functions #
    #################
    def density(self,savename,df,uv_name,optical_name):
        # Anna gave me a function for KDEs, I modified it to be a bit faster
        # Get rows where both columns exist
        rf = df[np.isfinite(df[uv_name]) & np.isfinite(df[optical_name])].reset_index(drop=True)
        # Get X and Y values
        x = (rf[uv_name] - rf[optical_name]).values
        y = rf[uv_name].values
        # Randomly sample to reduce computation time 
        random_index = np.random.choice(rf.shape[0],100000,replace=False)
        x = x[random_index]
        y = y[random_index]
        # Top cat styling for bulk of points 
        xy = np.vstack([x,y])
        z = stats.gaussian_kde(xy)(xy)
        index = z.argsort()
        x, y, z = np.array(x)[index], np.array(y)[index], np.array(z)[index]
        # Save data
        df.loc[random_index,f'{uv_name}_{optical_name}_x'] = x
        df.loc[random_index,f'{uv_name}_{optical_name}_y'] = y
        df.loc[random_index,f'{uv_name}_{optical_name}_z'] = z
        df.to_csv(savename,index=False)

    def truncate_colormap(self,cmap, minval=0.0, maxval=1.0, n=100):
        # Anna gave me this function to do top cat styling on kdes
        new_cmap = colors.LinearSegmentedColormap.from_list(
            'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
            cmap(np.linspace(minval, maxval, n)))
        return new_cmap
        
	
    def create_kde_files(self,data_dir):
        df = pd.read_csv(f'{data_dir}0_SUMS_Catalogs/CompleteCatalog/Crossmatched/{self.galaxy}_step5_crossmatch.csv') 

        # Drop if not a member
        df = df[df.gaia_pm_cut == "no"]
        df = df[df.gaia_px_cut == "no"].reset_index(drop=True)

        # Get rows where all columns exist 
        df = df[np.isfinite(df['uvm2_dered']) & np.isfinite(df['V_dered']) & np.isfinite(df['uvw2_dered']) & np.isfinite(df['U_dered']) & np.isfinite(df['uvw1_dered']) & np.isfinite(df['B_dered'])].reset_index(drop=True)

        # Only use columns we need 
        uv_cols = [f'{uv}_dered' for uv in ['uvm2','uvw2','uvw1']]
        optical_cols = [f'{opt}_dered' for opt in ['U','V','B']]
        color_cols = [f'{uv} - {opt} distance' for uv, opt in zip(['uvm2','uvw2','uvw1','u'],['v','u','b','v'])]
        
        # We also need to generate color columns for u - v  and uvw2 - u because we havent done that yet
        data_x = df['U_dered'] - df['V_dered']
        data_y = df['U_dered']
        #data_x_err = np.sqrt(df['e_U']**2 + df['e_V']**2)
        #color= self.AssessColor(data_x,data_y,data_x_err,self.zams['u'],self.zams['v'])
        df['u - v distance'] = self.DistanceToZams(data_x,data_y,self.zams['u'],self.zams['v'])
        
        data_x = df['uvw2_dered'] - df['U_dered']
        data_y = df['uvw2_dered']
        #data_x_err = np.sqrt(df['uvw2_err']**2 + df['e_U']**2)
       # color = self.AssessColor(data_x,data_y,data_x_err,self.zams['uvw2'],self.zams['u'])
        df['uvw2 - u distance'] = self.DistanceToZams(data_x,data_y,self.zams['uvw2'],self.zams['u'])
        
        # Rename for simplicity 
        df = df[uv_cols + optical_cols + color_cols]
        new_cols = ['uvm2','uvw2','uvw1','u','v','b','uvm2 - v distance','uvw2 - u distance','uvw1 - b distance','u - v distance']
        df.columns = new_cols 
        
        save_names = ['uvm2_v','uvw2_u','uvw1_b','u_v']
        for name in save_names:
            uv_name, optical_name = name.split('_')
            color_name = f'{uv_name} - {optical_name}'
            
            # Only get kde for things not very blue 
            ndf = df.copy()
            ndf = ndf[(ndf[f'{color_name} distance'] > 0)].reset_index(drop=True)

            # Randomly sample to reduce computation time 
            size = 100000
            if ndf.shape[0] > size:
                random_index = np.random.choice(ndf.shape[0],size,replace=False)
                ndf = ndf.iloc[random_index].reset_index(drop=True)

            # Get X and Y values
            x = ndf[uv_name].values - ndf[optical_name].values
            y = ndf[uv_name].values
        
            # Top cat styling for bulk of points 
            print(f'Running topcat style KDE for {name}')
            xy = np.vstack([x,y])
            z = stats.gaussian_kde(xy)(xy)
            index = z.argsort()
            x, y, z = np.array(x)[index], np.array(y)[index], np.array(z)[index]

            ndf[f'{name}_x'] = x
            ndf[f'{name}_y'] = y
            ndf[f'{name}_z'] = z

            ndf.to_csv(data_dir+f'0_SUMS_Catalogs/KDEs/{self.galaxy}_{name}_kde.csv',index=False)

    ###################
    # Bluewards of MS #
    ###################
    def AssessColor(self,data_x,data_y,data_x_err,zams_blue,zams_red):
        curve_x = np.array(zams_blue) - np.array(zams_red)
        curve_y = np.array(zams_blue)
        x_zams = np.interp(data_y,np.flip(curve_y,0),np.flip(curve_x,0))

        # Check if it is blue or not
        colors = []
        for i,x in enumerate(data_x):
            # If left of zams than consider it overlap 
            # unless we prove that it is blue within errors
            if x < (x_zams[i]):
                color = 'overlap'

            # If truly left of zams even within errors
            # than blue
            if x < (x_zams[i]-data_x_err[i]):
                color = 'blue'

            # If right of zams than consider it overlap
            # unless we prove that it is red within errors
            if x > x_zams[i]:
                color = 'overlap'

            # If right of zams than red
            if x > x_zams[i] + data_x_err[i]:
                color = 'red'

            # If color is off in narnia or a NaN
            if x < -90 or np.isnan(x):
                color = 'false'

            colors.append(color)
        return colors
    
    ####################
    # Distance From MS #
    ####################
    def DistanceToZams(self,data_x,data_y,zams_blue,zams_red):
        curve_x = np.array(zams_blue) - np.array(zams_red)
        curve_y = np.array(zams_blue)
        x_zams = np.interp(data_y,np.flip(curve_y,0),np.flip(curve_x,0))
        distance = data_x - x_zams
        # If distance is an extreme value make it nan 
        distance[np.abs(distance) > 20] = np.nan
        return distance


    def DistanceFromZamsForModels(self,df,blue_band,red_band):
        data_x = df[blue_band] - df[red_band]
        data_y = df[blue_band]
        zams_blue = self.zams[blue_band]
        zams_red = self.zams[red_band]
        curve_x = np.array(zams_blue) - np.array(zams_red)
        curve_y = np.array(zams_blue)
        f = interpolate.interp1d(curve_x, curve_y)
        xnew = np.linspace(np.min(curve_x),np.max(curve_x),1000)
        ynew = f(xnew)

        # Where is data_y closest to ynew?
        distances = []
        for ind, row in df.iterrows():
            x = row[blue_band] - row[red_band]
            y = row[blue_band]
            i_ymin = np.argmin(np.abs(ynew - y))
            i_x = xnew[i_ymin]
            d = x - i_x
            distances.append(d)

        return np.array(distances)

    ##############
    # Modify DFs # 
    ##############
    def format_candidates(self):
        self.candidates= self.candidates[self.candidates['galaxy'] == self.galaxy].reset_index(drop=True)
        # Generate color columns for u - v  and uvw2 - u 
        data_x = self.candidates['U_dered'] - self.candidates['V_dered']
        data_y = self.candidates['U_dered']
        data_x_err = np.sqrt(self.candidates['e_U']**2 + self.candidates['e_V']**2)
        
      
        color = self.AssessColor(data_x,data_y,data_x_err,self.zams['u'],self.zams['v'])
        self.candidates['u - v'] = color

        data_x = self.candidates['uvw2_dered'] - self.candidates['U_dered']
        data_y = self.candidates['uvw2_dered']
        data_x_err = np.sqrt(self.candidates['uvw2_err']**2 + self.candidates['e_U']**2)
        color = self.AssessColor(data_x,data_y,data_x_err,self.zams['uvw2'],self.zams['u'])
        self.candidates['uvw2 - u'] = color

        # Save just what we need ['uvm2_v','uvw2_u','uvw1_b','u_v']
        cols = ['ra','dec','galaxy','cut','discovery_name','uvw2_dered', 'uvw2_err', 'uvw1_dered', 'uvw1_err','uvm2_dered', 'uvm2_err', 'U_dered', 'e_U', 'B_dered', 'e_B','V_dered', 'e_V', 'I_dered', 'e_I',
                'uvw2 - b', 'uvw2 - v', 'uvw2 - i', 'uvw1 - b', 'uvw1 - v', 'uvw1 - i', 'uvm2 - b', 'uvm2 - v', 'uvm2 - i', 'u - v']
        self.candidates = self.candidates[cols]

        # Rename columns for simplicity
        self.candidates.columns = ['ra','dec','galaxy','cut','discovery_name','uvw2','uvw2_err','uvw1','uvw1_err','uvm2','uvm2_err','u','u_err','b','b_err','v','v_err','i','i_err',
                      'uvw2 - b', 'uvw2 - v', 'uvw2 - i', 'uvw1 - b', 'uvw1 - v', 'uvw1 - i', 'uvm2 - b', 'uvm2 - v', 'uvm2 - i', 'u - v']
                      
        
    def format_composites(self):
        photometry_columns = ['UVW2', 'UVM2', 'UVW1','U', 'B', 'V', 'I']
        # Rename just these columns
        for col in photometry_columns:
            self.composites.rename(columns={f'{col}':f'{col.lower()}'},inplace=True)

    ##############
    #  Plotting  # 
    ##############   

    def zams_threshold_plot(self, ax, blue_band, red_band, color = '#E0D4C4', overcolor = '#7E7875'):
        zams_x = self.zams[blue_band] - self.zams[red_band]
        zams_y = self.zams[blue_band]
        zams_y_lim = self.zams.loc[self.zams.m == 14.87, blue_band].values[0]
        
        # Interpolate or you get gaps in the zams line
        f = interpolate.interp1d(zams_x, zams_y)
        zams_x = np.linspace(np.min(zams_x), np.max(zams_x), 1000)
        zams_y = f(zams_x)
        diff = np.abs(zams_y - zams_y_lim)
        zams_y_lim = zams_y[diff == min(diff)][0]
        

        cmap = ListedColormap([color, overcolor])
        norm = BoundaryNorm([np.min(zams_y), zams_y_lim, np.max(zams_y)], cmap.N)
        points = np.array([zams_x, zams_y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)
        lc = LineCollection(segments, cmap=cmap, norm=norm, zorder=3)
        lc.set_array(zams_y)
        ax.add_collection(lc)
        ax.set_xlim(np.min(zams_x), np.max(zams_x))
        ax.set_ylim(np.min(zams_y) * 1.1, np.max(zams_y) * 1.1)
        return lc,zams_x,zams_y

    def plot_arrow(self,ax,x,y,blue_band,red_band,offset_x, offset_y, offset_r, fontsize,color='k'):
        # Deredden a point
        red_uv = y
        red_opt = y - x
        dered_uv = self.deredden(red_uv,blue_band)
        dered_opt = self.deredden(red_opt,red_band)
        # Get starting coordinate
        arrow_x = red_uv - red_opt
        arrow_y = red_uv
        # Get ending coordinate 
        arrow_dx = (dered_uv - dered_opt) - (red_uv - red_opt) 
        arrow_dy = dered_uv - red_uv
        # Hypotenuse
        hyp_x = arrow_x + arrow_dx
        hyp_y = arrow_y + arrow_dy
        # Adjacent
        adj_x = arrow_x
        adj_y = hyp_y
        # Right Triangle Lengths
        hypot_len = np.sqrt((hyp_x - arrow_x)**2 + (hyp_y - arrow_y)**2)
        adj_len = np.sqrt((adj_x - arrow_x)**2 + (adj_y - arrow_y)**2)
        # Get angle
        alpha = np.arccos(adj_len/hypot_len)
        # Plot arrow
        ax.arrow(arrow_x,arrow_y,arrow_dx,arrow_dy,head_width=0.1, head_length=0.1, fc=color, ec=color,zorder=10)
        ax.annotate(f'A$_v$ = {self.av}', xy=(arrow_dx/2 + arrow_x, arrow_dy/2 + arrow_y), 
                    xytext=(offset_x, offset_y), textcoords='offset points',
                    rotation= -90 + np.rad2deg(alpha) + offset_r,
                    color=color, family="serif",font="Times New Roman",fontsize=fontsize)
        
    def zams_box(self,ax,x,y,rotation,text,fontsize,textcolor,facecolor,edgecolor,zorder=20):
        props = dict(boxstyle='round', facecolor=facecolor, edgecolor=edgecolor, alpha=0.95)
        ax.text(x,y,text,fontsize=fontsize,zorder=zorder,weight='bold',rotation = rotation,color=textcolor, bbox=props,ha='center', va='center')

    def ostar_box(self,ax,x,rotation,y=14.3,text="O-type MS ",
                  fontsize=12,textcolor='dimgray',facecolor='#E0D4C4',edgecolor='silver',zorder=20):
        props = dict(boxstyle='round', facecolor=facecolor, edgecolor=edgecolor, alpha=0.95)
        ax.text(x,y,text,fontsize=fontsize,zorder=zorder,weight='bold',rotation = rotation,color=textcolor, bbox=props,ha='center', va='center',clip_on=True)

    def bstar_box(self,ax,x,rotation,y=17.5,text="B-type MS ",fontsize=12,textcolor='white',facecolor='#7E7875',edgecolor='#64605d',zorder=20):
        props = dict(boxstyle='round', facecolor=facecolor, edgecolor=edgecolor, alpha=0.95)
        ax.text(x,y,text,fontsize=fontsize,zorder=zorder,weight='bold',rotation = rotation,color=textcolor, bbox=props,ha='center', va='center')

#self.mist_fnames = ['0022000M','0024000M','0026000M','0030000M','0036000M','0040000M', '0046000M','0050000M','0054000M','0060000M','0066000M','0074000M']
    def Generate_MS_Table(self):
        ms_evol_tracks = []
        self.vega_ab_table = pd.read_csv('https://waps.cfa.harvard.edu/MIST/BC_tables/zeropoints.txt',delimiter='\\s+')
        for mist_file in self.mist_fnames:
            # Get the associated UBVI filename
            track_file = self.mist_dir +'MIST/Rotation/'+ mist_file + '.track.eep'
            swift_file = self.mist_dir +'Swift/Rotation/'+ mist_file + '.track.eep.cmd'
            ubvi_file = self.mist_dir +'UBVI/Rotation/'+ mist_file + '.track.eep.cmd'
            # Get MIST data
            mist = self.Open_File(track_file,11,mist=True)
            swift = self.Open_File(swift_file,14,mist=True)
            ubvi = self.Open_File(ubvi_file,14,mist=True)
            # Get the rows  after some frac_MS age
            start_age = self.Get_Age(mist,frac_MS=0.)
            stop_age = self.Get_Age(mist,frac_MS=1.1)
            mist = self.Get_Row_By_Age(mist,start_age,stop_age)
            swift = self.Get_Row_By_Age(swift,start_age,stop_age)
            ubvi = self.Get_Row_By_Age(ubvi,start_age,stop_age)

            # Get the rows
            rows = {
                'age':np.round(mist.star_age,2),
                # Swift is already in AB
                'UVW2':self.Absolute_to_Apparent(swift.Swift_UVW2,self.distance),
                'UVM2':self.Absolute_to_Apparent(swift.Swift_UVM2,self.distance),
                'UVW1':self.Absolute_to_Apparent(swift.Swift_UVW1,self.distance),
                # Convert from Vega to AB for UBVI
                'U':self.Absolute_to_Apparent(self.Vega2AB(ubvi,'Bessell_U'),self.distance),
                'B':self.Absolute_to_Apparent(self.Vega2AB(ubvi,'Bessell_B'),self.distance),
                'V':self.Absolute_to_Apparent(self.Vega2AB(ubvi,'Bessell_V'),self.distance),
                'I':self.Absolute_to_Apparent(self.Vega2AB(ubvi,'Bessell_I'),self.distance)}

            ms_evol_tracks.append(pd.DataFrame(rows).reset_index(drop=True))

        return ms_evol_tracks


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
        return age
    
    def Get_Row_By_Age(self,df,start_age,stop_age):
        # Find where frac_MS = 0.2
        age_diff = np.abs(df['star_age'] - start_age)
        start_ind = np.argmin(age_diff)
        age_diff = np.abs(df['star_age'] - stop_age)
        stop_ind = np.argmin(age_diff)

        rows =  df.iloc[start_ind:stop_ind]
        return rows
    
    def Vega2AB(self,df,band):
        # Swift is already in AB but MIST needs to be converted from Vega
        vega_ab = self.vega_ab_table.loc[self.vega_ab_table['filter'] == band,'mag(Vega/AB)']
        converted_mag = df[band] + vega_ab.values[0]
        return converted_mag
    
    def Absolute_to_Apparent(self,AbsoluteMag,distance):
        # Convert Absolute to Apparent Magnitude
        return AbsoluteMag + 5 * (np.log10(distance/10))
    