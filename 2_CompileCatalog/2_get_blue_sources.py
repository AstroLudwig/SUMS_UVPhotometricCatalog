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
galaxy = 'lmc'
path = f'Step4/{galaxy}_photometry.csv'
save = f'Step5/{galaxy}_colors.csv'
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

#Function to figure out bluewards vs. not: 
# Ex: Zams_blue = 'uvm2', zams_red = 'v' for colors 'uvm2 - v'
def AssessColor(data_x,data_y,data_x_err,zams_blue,zams_red):

	curve_x = np.array(zams_blue) - np.array(zams_red)
	curve_y = np.array(zams_blue)
	x_zams = np.interp(data_y,np.flip(curve_y,0),np.flip(curve_x,0))
	#ZAMS_distance = data_x - zams_color_ref

	#Have a boolean flag for if it is close. 
	colors,overlap = [],[]
	for i,x in enumerate(data_x):
		#print(x,x_zams[i])
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
# For SMC 
if galaxy == 'smc':
	distance = smc_distance
	# Goetberg and MIST Models
	models = {
		"zams_vega": pd.read_csv('../1_run_tractor/models/ZAMS_Z0.002_Vegamag.txt',sep='\s+',comment='#'),
		"shs_vega":  pd.read_csv('../1_run_tractor/models/stripped_stars_Z0.002_Vegamag.txt',sep='\s+',comment='#'),
		"zams_ab":   pd.read_csv('../1_run_tractor/models/ZAMS_Z0.002_ABmag.txt',sep='\s+',comment='#'),
		"shs_ab":	 pd.read_csv('../1_run_tractor/models/stripped_stars_Z0.002_ABmag.txt',sep='\s+',comment='#')		
	}
	# Extinction correctionscalculated here: /FoesHackSession/Hack1 <---- ask maria about this
	extinction_coeff = {
		'uvw2_mag' : 3.522767303652845,
		'uvm2_mag' : 3.062652780996107,
		'uvw1_mag' : 2.7436011426496876,
		'U'  : 1.746851127566682,
		'B'  : 1.4073996800012445,
		'V'  : 1.0353852912271932,
		'I'  : 0.5911705440475679,
	}
# For LMC
if galaxy == 'lmc':
	distance = lmc_distance
	# Goetberg and MIST Models
	models = {
		"zams_vega": pd.read_csv('../1_run_tractor/models/ZAMS_Z0.006_Vegamag.txt',sep='\s+',comment='#'),
		"shs_vega":  pd.read_csv('../1_run_tractor/models/stripped_stars_Z0.006_Vegamag.txt',sep='\s+',comment='#'),
		"zams_ab":   pd.read_csv('../1_run_tractor/models/ZAMS_Z0.006_ABmag.txt',sep='\s+',comment='#'),
		"shs_ab":	 pd.read_csv('../1_run_tractor/models/stripped_stars_Z0.006_ABmag.txt',sep='\s+',comment='#')
	}
	# Extinction correctionscalculated here: /FoesHackSession/Hack1 <---- ask maria about this
	extinction_coeff = {
		'uvw2_mag' : 2.644541680555541,
		'uvm2_mag' : 2.7233599880955124,
		'uvw1_mag' : 2.3536449902431045,
		'U'  : 1.5634790597438197,
		'B'  : 1.3170082045312625,
		'V'  : 1.0333402844940784,
		'I'  : 0.7366865234305091,
	}
# Generalized / For both lmc and smc #
# MIST zams model in apparent mags
zams = {
	'uvw2':   AbsoluteToApparent(models['zams_ab']["UVW2_spec"],distance),
	'uvw1':   AbsoluteToApparent(models['zams_ab']["UVW1_spec"],distance),
	'uvm2':   AbsoluteToApparent(models['zams_ab']["UVM2_spec"],distance),
	'u':      AbsoluteToApparent(models['zams_ab']["U_spec"],distance),
	'b':      AbsoluteToApparent(models['zams_ab']["B_spec"],distance),
	'v':      AbsoluteToApparent(models['zams_ab']["V_spec"],distance),
	'i':         AbsoluteToApparent(models['zams_ab']["I_spec"],distance)
}
# Stripped helium star models in apparent mags
shs = {  
	'uvm2_ab': AbsoluteToApparent(models['shs_ab']["UVM2"],distance),
	'v_ab':    AbsoluteToApparent(models['shs_ab']["V"],distance),
}
# Conversion from vega to ab coefficient used in Gotberg's models. 
vega2ab = { 
	'uvw2_mag' : (models['shs_ab']['UVW2'] - models['shs_vega']['UVW2'])[0],
	'uvm2_mag' : (models['shs_ab']['UVM2'] - models['shs_vega']['UVM2'])[0],
	'uvw1_mag' : (models['shs_ab']['UVW1'] - models['shs_vega']['UVW1'])[0],
	'U' :    (models['shs_ab']['U'] - models['shs_vega']['U'])[0],
	'B' :    (models['shs_ab']['B'] - models['shs_vega']['B'])[0],
	'V' :    (models['shs_ab']['V'] - models['shs_vega']['V'])[0],
	'I' :    (models['shs_ab']['I'] - models['shs_vega']['I'])[0],
	'J' :    (models['shs_ab']['J'] - models['shs_vega']['J'])[0],
	'H' :    (models['shs_ab']['H'] - models['shs_vega']['H'])[0],
	'Ks' :    (models['shs_ab']['Ks'] - models['shs_vega']['Ks'])[0],
}
print(vega2ab)
print(f'{galaxy} models loaded')
###############
#  Photometry #
###############   

df = pd.read_csv(path)

# Replace null optical with NaN 
for uvfilter,error in zip(['U','B','V','I'],['e_U','e_B','e_V','e_I']):
	# For no extincted photometry
	df[f'{uvfilter[:4]}'].replace(0.0,np.nan)
	# For no photometric errors
	df.loc[df[error] == 0,f'{uvfilter}'] = np.nan
    
print(f'{galaxy} photometry loaded')
####################################
# Combined errors on UV Magnitudes # - We don't end up using this 
####################################

# If std is positive then combine with photometry error in quadrature.
# If not just use the magnitude error. 

# for uvfilter in ['uvw2','uvw1','uvm2']:			 # Conditions
# 	# By Residual Frac
#     df[f'{uvfilter}_combined_err'] = np.select([(df[f'{uvfilter}_mag_std']>0),(df[f'{uvfilter}_mag_std']<0),(np.isnan(df[f'{uvfilter}_mag']))],
# 												 # Values
# 												[Quad(df[f'{uvfilter}_mag_err'],df[f'{uvfilter}_mag_std']),
#                                                  df[f'{uvfilter}_mag_err'],
#                                                  np.repeat(np.nan,len(np.isnan(df[f'{uvfilter}_mag_std'])))])

    
#print(f'{galaxy} combined errors')
####################################################
# Convert magnitudes from vega (via heasarc) to AB #  
#################################################### 

for uvfilter in ['uvw2_mag','uvw1_mag','uvm2_mag']:	

	df[uvfilter] = df[uvfilter] + vega2ab[uvfilter]

for mcpsfilter in ['U','B','V','I','J','H','Ks']:	

	df[mcpsfilter] = df[mcpsfilter] + vega2ab[mcpsfilter]
    

#####################################
# Correct magnitudes for extinction #
#####################################

if galaxy == 'smc':
	Av = 0.22
if galaxy == 'lmc':
	Av = 0.38

for uvfilter in ['uvw2_mag','uvm2_mag','uvw1_mag','U','B','V','I']:	
    
        df[f'{uvfilter}_dered'] = df[uvfilter] - Av * extinction_coeff[uvfilter]
        
print(f'{galaxy} extinction corrected')
################################
# Calculate colors and errors: #
################################


colors = {
	'uvw2 - b' : df['uvw2_mag_dered'] - df['B_dered'],
	'uvw2 - v' : df['uvw2_mag_dered'] - df['V_dered'],
	'uvw2 - i' : df['uvw2_mag_dered'] - df['I_dered'],
	'uvw1 - b' : df['uvw1_mag_dered'] - df['B_dered'],
	'uvw1 - v' : df['uvw1_mag_dered'] - df['V_dered'],
	'uvw1 - i' : df['uvw1_mag_dered'] - df['I_dered'],
	'uvm2 - b' : df['uvm2_mag_dered'] - df['B_dered'],
	'uvm2 - v' : df['uvm2_mag_dered'] - df['V_dered'],
	'uvm2 - i' : df['uvm2_mag_dered'] - df['I_dered']
}

color_errs = {
	'uvw2 - b' : Quad(df['uvw2_mag_err'] , df['e_B']),
	'uvw2 - v' : Quad(df['uvw2_mag_err'] , df['e_V']),
	'uvw2 - i' : Quad(df['uvw2_mag_err'] , df['e_I']),
	'uvw1 - b' : Quad(df['uvw1_mag_err'] , df['e_B']),
	'uvw1 - v' : Quad(df['uvw1_mag_err'] , df['e_V']),
	'uvw1 - i' : Quad(df['uvw1_mag_err'] , df['e_I']),
	'uvm2 - b' : Quad(df['uvm2_mag_err'] , df['e_B']),
	'uvm2 - v' : Quad(df['uvm2_mag_err'] , df['e_V']),
	'uvm2 - i' : Quad(df['uvm2_mag_err'] , df['e_I'])
}



##############################
# Calculate if blue of zams: #
##############################

color_labels = {                   # x             # y                   # x_err       # zams_blue  # zams_red
	'uvw2 - b' : AssessColor(colors['uvw2 - b'],df['uvw2_mag_dered'],color_errs['uvw2 - b'],zams['uvw2'],zams['b'])[0],
	'uvw2 - v' : AssessColor(colors['uvw2 - v'],df['uvw2_mag_dered'],color_errs['uvw2 - v'],zams['uvw2'],zams['v'])[0],
	'uvw2 - i' : AssessColor(colors['uvw2 - i'],df['uvw2_mag_dered'],color_errs['uvw2 - i'],zams['uvw2'],zams['i'])[0],
	'uvw1 - b' : AssessColor(colors['uvw1 - b'],df['uvw1_mag_dered'],color_errs['uvw1 - b'],zams['uvw1'],zams['b'])[0],
	'uvw1 - v' : AssessColor(colors['uvw1 - v'],df['uvw1_mag_dered'],color_errs['uvw1 - v'],zams['uvw1'],zams['v'])[0],
	'uvw1 - i' : AssessColor(colors['uvw1 - i'],df['uvw1_mag_dered'],color_errs['uvw1 - i'],zams['uvw1'],zams['i'])[0],
	'uvm2 - b' : AssessColor(colors['uvm2 - b'],df['uvm2_mag_dered'],color_errs['uvm2 - b'],zams['uvm2'],zams['b'])[0],
	'uvm2 - v' : AssessColor(colors['uvm2 - v'],df['uvm2_mag_dered'],color_errs['uvm2 - v'],zams['uvm2'],zams['v'])[0],
	'uvm2 - i' : AssessColor(colors['uvm2 - i'],df['uvm2_mag_dered'],color_errs['uvm2 - i'],zams['uvm2'],zams['i'])[0],

	'uvw2 - b overlap' : AssessColor(colors['uvw2 - b'],df['uvw2_mag_dered'],color_errs['uvw2 - b'],zams['uvw2'],zams['b'])[1],
	'uvw2 - v overlap' : AssessColor(colors['uvw2 - v'],df['uvw2_mag_dered'],color_errs['uvw2 - v'],zams['uvw2'],zams['v'])[1],
	'uvw2 - i overlap' : AssessColor(colors['uvw2 - i'],df['uvw2_mag_dered'],color_errs['uvw2 - i'],zams['uvw2'],zams['i'])[1],
	'uvw1 - b overlap' : AssessColor(colors['uvw1 - b'],df['uvw1_mag_dered'],color_errs['uvw1 - b'],zams['uvw1'],zams['b'])[1],
	'uvw1 - v overlap' : AssessColor(colors['uvw1 - v'],df['uvw1_mag_dered'],color_errs['uvw1 - v'],zams['uvw1'],zams['v'])[1],
	'uvw1 - i overlap' : AssessColor(colors['uvw1 - i'],df['uvw1_mag_dered'],color_errs['uvw1 - i'],zams['uvw1'],zams['i'])[1],
	'uvm2 - b overlap' : AssessColor(colors['uvm2 - b'],df['uvm2_mag_dered'],color_errs['uvm2 - b'],zams['uvm2'],zams['b'])[1],
	'uvm2 - v overlap' : AssessColor(colors['uvm2 - v'],df['uvm2_mag_dered'],color_errs['uvm2 - v'],zams['uvm2'],zams['v'])[1],
	'uvm2 - i overlap' : AssessColor(colors['uvm2 - i'],df['uvm2_mag_dered'],color_errs['uvm2 - i'],zams['uvm2'],zams['i'])[1]

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