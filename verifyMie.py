import os
import argparse
import numpy as np
import miepython
from visualize import plt_qcoeffs
from verifyEpsilon import permittivity_ccpr
from constants import FILE_QFDTD_MIE, FILE_SIM, FILE_CCPR
import IO 
import sys

'''
Comparison of the absorption and scattering efficiencies of the Mie theory with the results of the 
in-house FDTD software with TopOpt implementation. The program analyze.py should
be run first for a Sphere to obtain the Qabs values obtained from the volume integration.
Mie's cross sections are computed using the library 'miepython',
-> https://miepython.readthedocs.io/en/latest/02_efficiencies.html
'''

def verify_mie(directory):
	
	# Get radius which has been used in FDTD sim
	file_path = os.path.join(directory, FILE_SIM)
	radius = IO.readin_radius(file_path)		
	print("Mie Sphere radius [m]:", "{:.1e}".format(radius))
	# Get wavelenght range based on FDTD code input file
	lambd = IO.readin_wavelength(file_path)

	# Read in CCPR parameters and compute permittivity
	file_path = os.path.join(directory, FILE_CCPR)
	p_ccpr = IO.readin_CCPR_params(file_path)
	lambd, eps_ccpr = permittivity_ccpr(p_ccpr, lambd, True)
	
	# Verify by visualizing the permittivity
	# visualize.plt_refractive_index(lambd, eps_ccpr, timeconvention=-1)
	n_ccpr = np.sqrt(eps_ccpr)

	isotropic = np.all(n_ccpr == n_ccpr[0, :], axis=0).all()
	if not isotropic:
		raise ValueError("Permittivity/refractive index is anisotropic. Mie Theory comparison is not covered for that yet!")
	else:
		eps_ccpr = eps_ccpr[0,:]
		n_ccpr = n_ccpr[0,:]

	# Compute coeffients based on Mie's theory
	x = 2 * np.pi * radius / lambd
	qext, qsca, qback, g = miepython.mie(n_ccpr, x)
	qabs = qext - qsca
	'''
	Readin coefficients (TFSF surface integr.) based on FDTD code output
	After transposing table:
		dim=0 : wavelenght (in units nm!)
		dim=1 : Qscatt
		dim=2 : Qabs
		dim=3 : Qext
	'''
	# Flip transposed array to reorder quantities from smalles to largest wavelength
	lambd_data, qsca_data, qabs_data, qext_data = np.flip(np.transpose(np.loadtxt(os.path.join(directory, FILE_QFDTD_MIE))), axis=1)
	# Convert wavelenght in nm to m
	lambd_data = 1e-9 * lambd_data 

	# Check if lambda range is same for Mie and data (for comparison)
	print("Wl start: \nMie: ", "{:.1e}".format(lambd[0]),  "\nData:", "{:.1e}".format(lambd_data[0])) 
	print("Wl end:   \nMie: ", "{:.1e}".format(lambd[-1]), "\nData:", "{:.1e}".format(lambd_data[-1]))
	
	# Replace the Qabs computed by surface integral (in in-house Code) method with the volume integral method
	file_str_abs_integr = os.path.join(directory, "lamb_qabs_radius.npy")
	try:
		lambd_integr, qabs_integr, radius_integr = np.load(file_str_abs_integr, allow_pickle=True)
	except:
		print("\n(!) file ", file_str_abs_integr, " is missing!") 
		print("In order to display a comparison of Qabs \n" \
		"of a sphere between Mie theory and the FDTD simulation (with volume integr.), \n" \
		"the program analyze.py must first be executed for a sphere.") 
		print("Program is terminated, sorry")
		sys.exit()	

	plt_qcoeffs(os.path.join(directory), lambd, qext, qsca, qabs, lambd_data, qext_data, qsca_data, qabs_data, lambd_integr, qabs_integr)

	return radius, lambd, qabs, lambd_data, qabs_data
	
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('--dir_mie', required=True, type=str, help='Directory from which the data is read in')
	args = parser.parse_args()
	verify_mie(args.dir_mie)