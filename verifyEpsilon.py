import os
import numpy as np
import argparse
import IO
import visualize
from constants import V_LIGHT, EPS_0

SKIPROWS = 10 # Rows skipped in exp. data file

def permittivity_ccpr(p, lambd, endpoint=True):
    """
    Returns wavelength array and (component-wise) permittivity based on the CCPR model 
    within a given wavelenght range.

    Parameters
    ----------
    p: dictionary 
        Dictionary containing CCPR parameters.
    lambd: array 
        Array containing wavelength in units [m].
    """
    omega = 2 * np.pi * V_LIGHT / lambd

    poles = p['c'][:,:,np.newaxis] / (1j*omega-p['a'][:,:,np.newaxis]) + \
            np.conjugate(p['c'][:,:,np.newaxis]) / (1j*omega - np.conjugate(p['a'][:,:,np.newaxis]))
    eps_CCPR = p['eps'][:, np.newaxis] + p['sigm'][:, np.newaxis] / (1j * omega * EPS_0) + np.sum(poles, axis=1) 

    return lambd, eps_CCPR

def verify_epsilon(lambd_micron, dir_CCPR, dir_data):
    """
    Comparison of the CCPR model fit with experimental data.
    """
    # ----------------------------------------------------------
    # Compute epsilon based on CCPR parameters
    # ----------------------------------------------------------
    p_ccpr = IO.readin_CCPR_params(dir_CCPR)
    lambd = 1e-6 * lambd_micron # microns -> m
    lambd_ccpr, eps_ccpr = permittivity_ccpr(p_ccpr, lambd, True)
    # ----------------------------------------------------------
    # Readin the experimental data
    # ----------------------------------------------------------
    data = np.loadtxt(dir_data, skiprows=SKIPROWS).flatten()
    data = np.reshape(data, (int(len(data) / 3), 3))
    lambd_data_micron, n, k  = data[:,0], data[:,1], data[:,2]
    idxUp = np.where((lambd_data_micron > lambd_micron[0]))[0]
    idxLow = np.where((lambd_data_micron < lambd_micron[-1]))[0]
    idxLimits = np.intersect1d(idxUp , idxLow)
    eps_data = (n + 1j * k)[idxLimits]**2
    lambd_data_m_bounded = 1e-6 * lambd_data_micron[idxLimits]
    lambd_data_m = 1e-6 * lambd_data_micron[idxLimits]
    # ----------------------------------------------------------
    # Visualize
    # ----------------------------------------------------------
    visualize.plt_refractive_index(os.path.dirname(dir_CCPR), lambd_ccpr, eps_ccpr, -1, True, lambd_data_m, eps_data, True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--lower_bound', type=float, required=True, help="The lower wavelength boundary value in microns.")
    parser.add_argument('--upper_bound', type=float, required=True, help="The upper wavelength boundary value in microns.")
    parser.add_argument('--dir_CCPR', type=str, required=True, help='Directory from which the CCPR parameters are read in')
    parser.add_argument('--dir_data', type=str, help='Directory from which the experimental data is read in')
    parser.add_argument('--N', required=True, type=int, help='Number of wavelength points')
    args = parser.parse_args()
    lambd_limits = np.array([args.lower_bound, args.upper_bound]) # in microns
    lambd = np.linspace(lambd_limits[0], lambd_limits[-1], args.N, endpoint=True)
    verify_epsilon(lambd, args.dir_CCPR, args.dir_data)