import sys 
import glob
sys.path.insert(0, 'lib/')
from Heasarc import *


def run_heasarc(path):
	#########################################
	def format_filter(uvfilter):
		if uvfilter == 'um2':
			return 'uvm2'
		if uvfilter == 'uw2':
			return 'uvw2'
		if uvfilter == 'uw1':
			return 'uvw1'
		else:
			print('Appropriate filter not found.')
			return None


	n = len(path)
	folders = glob.glob(f'{path}*')
	print(f"{len(folders)} files to process.")
	for f in folders:
		obsid = f[n+5:n+10]
		uvfilter = f[n+13:n+16]
		uv_filter = format_filter(uvfilter)
		segment = f[-3]
		extension = f[-1]
		tractor_file = f + f'/{obsid}_{uvfilter}_{segment}_{extension}.fits'
		print(f"Running Heasarc on Obsid: {obsid} UVfilter:{uvfilter}")
		try:
			Magnitudes = HeasarcRoutines(tractor_file,uv_filter)
		except:
			print(f"Issue with {f}")
			continue
	print("Run Complete")

#path = f'/home/bethany/Desktop/SUMS_Tractor_Data/lmc/lmc_4558X/' # Change Path
path = '/home/bethany/Desktop/SUMS_Tractor_Data/lmc_missing_file/'
run_heasarc(path)
