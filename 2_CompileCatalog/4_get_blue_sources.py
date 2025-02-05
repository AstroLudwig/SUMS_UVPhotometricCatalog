import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np 
from scipy.integrate import simps
import time
# Create a version of the catalog with extinction and AB mag accounted for.
# Corrects the photometry for extinction and switches it into AB magnitudes. 
# At that point it also goes through and figures out if the photometry is 
# bluewards of the ZAMS in nine different UV optical CMDs.

# Choose which galaxy to work on: 
galaxy = "smc"
data_dir = "/home/bethany/Projects/0_Data"
step4_dir = data_dir + "0_SUMS_Catalogs/CompleteCatalog/Step4/"
step5_dir = data_dir + "0_SUMS_Catalogs/CompleteCatalog/Step5/"
path = step4_dir + f'{galaxy}_photometry.csv'
save = step5_dir + f'{galaxy}_colors.csv'

##############################
# Constants and definitions: #
##############################

smc_distance = 60.6e3
lmc_distance = 50e3

#############
# Functions #
#############
# Convert absolute mags to apparent
def AbsoluteToApparent(AbsoluteMag,distance):
	return AbsoluteMag + 5 * (np.log10(distance/10))

# Summing in quadrature
def Quad(X,Y):
	return np.sqrt(X**2+Y**2) 

# Function to figure out bluewards vs. not: 
# Ex: Zams_blue = 'uvm2', zams_red = 'v' for colors 'uvm2 - v'
def AssessColor(data_x,data_y,data_x_err,zams_blue,zams_red):
	curve_x = np.array(zams_blue) - np.array(zams_red)
	curve_y = np.array(zams_blue)
	x_zams = np.interp(data_y,np.flip(curve_y,0),np.flip(curve_x,0))

	# Check if it is blue or not
	colors,overlap = [],[]
	for i,x in enumerate(data_x):
		# If left of zams than blue
		if x < x_zams[i]:
			color = 'blue'

		# If right of zams than red
		elif x > x_zams[i]:
			color = 'red'

		# If color is off in narnia 
		elif x < -90:
			color = 'false'

		# If color is a NAN
		elif np.isnan(x):
			color = 'false'
			
		# Should anything else happen, don't think it will
		else:
			color = 'oops'

		colors.append(color)

		# If near zams within error than overlap 
		if (x > x_zams[i] - data_x_err[i]) and (x < x_zams[i] + data_x_err[i]):
			overlap.append('yes')
		else:
			overlap.append('no')
	return colors,overlap

###############
# M o d e l s #
###############
start = time.time()
# SMC Models 
if galaxy == 'smc':
	distance = smc_distance
	# Goetberg and MIST Models
	models = {
		"zams_ab":   pd.read_csv('C:/Projects/0_Data/Models/ZAMS/ZAMS_Z0.002_ABmag.txt',sep=r'\s+',comment='#'),
		"shs_ab":	 pd.read_csv('C:/Projects/0_Data/Models/SHS/stripped_stars_Z0.002_ABmag.txt',sep=r'\s+',comment='#')		
	}

# LMC Models
if galaxy == 'lmc':
	distance = lmc_distance
	# Goetberg and MIST Models
	models = {
		"zams_ab":   pd.read_csv('C:/Projects/0_Data/Models/ZAMS/ZAMS_Z0.006_ABmag.txt',sep=r'\s+',comment='#'),
		"shs_ab":	 pd.read_csv('C:/Projects/0_Data/Models/SHS/stripped_stars_Z0.006_ABmag.txt',sep=r'\s+',comment='#')
	}

# Convert to apparent mags
zams = {
	'uvw2':   AbsoluteToApparent(models['zams_ab']["UVW2_spec"],distance),
	'uvw1':   AbsoluteToApparent(models['zams_ab']["UVW1_spec"],distance),
	'uvm2':   AbsoluteToApparent(models['zams_ab']["UVM2_spec"],distance),
	'u':      AbsoluteToApparent(models['zams_ab']["U_spec"],distance),
	'b':      AbsoluteToApparent(models['zams_ab']["B_spec"],distance),
	'v':      AbsoluteToApparent(models['zams_ab']["V_spec"],distance),
	'i':         AbsoluteToApparent(models['zams_ab']["I_spec"],distance)
}
shs = {  
	'uvm2_ab': AbsoluteToApparent(models['shs_ab']["UVM2"],distance),
	'v_ab':    AbsoluteToApparent(models['shs_ab']["V"],distance),
}
print(f'{galaxy} models loaded')

###############
#  Photometry #
###############   

df = pd.read_csv(path)
# Remove _mag
df.columns = [col.replace('_mag','') for col in df.columns]

# Extinction and Vega to AB
ref = pd.read_csv('C:/Projects/0_Data/SUMS_Catalogs/ExtinctionAndVega2AB.csv')

# Replace null optical with NaN 
for uvfilter,error in zip(['U','B','V','I'],['e_U','e_B','e_V','e_I']):
	# For no extincted photometry
	df[f'{uvfilter[:4]}'].replace(0.0,np.nan)
	# For no photometric errors
	df.loc[df[error] == 0,f'{uvfilter}'] = np.nan
    
print(f'{galaxy} photometry loaded')

######################################
# Convert magnitudes from Vega to AB #  
######################################

for uvfilter in ['uvw2','uvm2','uvw1','U','B','V','I']:	

	df[uvfilter] = df[uvfilter] + ref.loc[ref['filter'] == uvfilter.upper(),'Vega to AB'].values[0]

    
#####################################
# Correct magnitudes for extinction #
#####################################

for uvfilter in ['uvw2','uvm2','uvw1','U','B','V','I']:	
    
        df[f'{uvfilter}_dered'] = df[uvfilter] - ref.loc[ref['filter'] == uvfilter.upper(),f'A_lambda ({galaxy.upper()})'].values[0]
        
print(f'{galaxy} extinction corrected')

################################
# Calculate colors and errors: #
################################
colors = {
	'uvw2 - b' : df['uvw2_dered'] - df['B_dered'],
	'uvw2 - v' : df['uvw2_dered'] - df['V_dered'],
	'uvw2 - i' : df['uvw2_dered'] - df['I_dered'],
	'uvw1 - b' : df['uvw1_dered'] - df['B_dered'],
	'uvw1 - v' : df['uvw1_dered'] - df['V_dered'],
	'uvw1 - i' : df['uvw1_dered'] - df['I_dered'],
	'uvm2 - b' : df['uvm2_dered'] - df['B_dered'],
	'uvm2 - v' : df['uvm2_dered'] - df['V_dered'],
	'uvm2 - i' : df['uvm2_dered'] - df['I_dered']
}

color_errs = {
	'uvw2 - b' : Quad(df['uvw2_err'] , df['e_B']),
	'uvw2 - v' : Quad(df['uvw2_err'] , df['e_V']),
	'uvw2 - i' : Quad(df['uvw2_err'] , df['e_I']),
	'uvw1 - b' : Quad(df['uvw1_err'] , df['e_B']),
	'uvw1 - v' : Quad(df['uvw1_err'] , df['e_V']),
	'uvw1 - i' : Quad(df['uvw1_err'] , df['e_I']),
	'uvm2 - b' : Quad(df['uvm2_err'] , df['e_B']),
	'uvm2 - v' : Quad(df['uvm2_err'] , df['e_V']),
	'uvm2 - i' : Quad(df['uvm2_err'] , df['e_I'])
}


##############################
# Calculate if blue of zams: #
##############################

color_labels = {                   # x             # y                   # x_err       # zams_blue  # zams_red
	'uvw2 - b' : AssessColor(colors['uvw2 - b'],df['uvw2_dered'],color_errs['uvw2 - b'],zams['uvw2'],zams['b'])[0],
	'uvw2 - v' : AssessColor(colors['uvw2 - v'],df['uvw2_dered'],color_errs['uvw2 - v'],zams['uvw2'],zams['v'])[0],
	'uvw2 - i' : AssessColor(colors['uvw2 - i'],df['uvw2_dered'],color_errs['uvw2 - i'],zams['uvw2'],zams['i'])[0],
	'uvw1 - b' : AssessColor(colors['uvw1 - b'],df['uvw1_dered'],color_errs['uvw1 - b'],zams['uvw1'],zams['b'])[0],
	'uvw1 - v' : AssessColor(colors['uvw1 - v'],df['uvw1_dered'],color_errs['uvw1 - v'],zams['uvw1'],zams['v'])[0],
	'uvw1 - i' : AssessColor(colors['uvw1 - i'],df['uvw1_dered'],color_errs['uvw1 - i'],zams['uvw1'],zams['i'])[0],
	'uvm2 - b' : AssessColor(colors['uvm2 - b'],df['uvm2_dered'],color_errs['uvm2 - b'],zams['uvm2'],zams['b'])[0],
	'uvm2 - v' : AssessColor(colors['uvm2 - v'],df['uvm2_dered'],color_errs['uvm2 - v'],zams['uvm2'],zams['v'])[0],
	'uvm2 - i' : AssessColor(colors['uvm2 - i'],df['uvm2_dered'],color_errs['uvm2 - i'],zams['uvm2'],zams['i'])[0],

	'uvw2 - b overlap' : AssessColor(colors['uvw2 - b'],df['uvw2_dered'],color_errs['uvw2 - b'],zams['uvw2'],zams['b'])[1],
	'uvw2 - v overlap' : AssessColor(colors['uvw2 - v'],df['uvw2_dered'],color_errs['uvw2 - v'],zams['uvw2'],zams['v'])[1],
	'uvw2 - i overlap' : AssessColor(colors['uvw2 - i'],df['uvw2_dered'],color_errs['uvw2 - i'],zams['uvw2'],zams['i'])[1],
	'uvw1 - b overlap' : AssessColor(colors['uvw1 - b'],df['uvw1_dered'],color_errs['uvw1 - b'],zams['uvw1'],zams['b'])[1],
	'uvw1 - v overlap' : AssessColor(colors['uvw1 - v'],df['uvw1_dered'],color_errs['uvw1 - v'],zams['uvw1'],zams['v'])[1],
	'uvw1 - i overlap' : AssessColor(colors['uvw1 - i'],df['uvw1_dered'],color_errs['uvw1 - i'],zams['uvw1'],zams['i'])[1],
	'uvm2 - b overlap' : AssessColor(colors['uvm2 - b'],df['uvm2_dered'],color_errs['uvm2 - b'],zams['uvm2'],zams['b'])[1],
	'uvm2 - v overlap' : AssessColor(colors['uvm2 - v'],df['uvm2_dered'],color_errs['uvm2 - v'],zams['uvm2'],zams['v'])[1],
	'uvm2 - i overlap' : AssessColor(colors['uvm2 - i'],df['uvm2_dered'],color_errs['uvm2 - i'],zams['uvm2'],zams['i'])[1]

}
print(f'{galaxy} colors calculated')

########################
# Append dictionaries: #
########################

for key in color_labels.keys():

	df[key] = color_labels[key]

#########
# Save: #
#########
print(f'{galaxy} saving...')
df.to_csv(save)

print(f'{galaxy} complete. Time taken: {time.time() - start} s')