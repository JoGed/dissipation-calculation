import argparse
import os
import IO
import visualize
import numpy as np
from matplotlib import pyplot as plt
from constants import V_LIGHT, EPS_0, FILE_CCPR, FILE_SIM
import re

"""
    Script that performs the computation of the electrical power dissipation 
    density based on the electric field data and the complex permittivity.
"""

def topt_regions_idx(p_sim):
    """
    Returns the indices of design and obervation region within
    the index scheme used for the original DFT fields. 
    """
    dr_bounds = p_sim['dr'] - np.expand_dims(p_sim['cpml'], -1)
    or_bounds = p_sim['or'] - np.expand_dims(p_sim['cpml'], -1)
    or_in_dr_bounds = or_bounds - np.expand_dims(dr_bounds[:, 0], -1)

    return dr_bounds, or_bounds, or_in_dr_bounds

def interpolate_components(arr, ftype="vector"):
    """
    Returns an interpolation of the scalar/field data originally defined
    on the staggered Yee grid in 3D.

    Parameters
    ----------
    arr : array
        Values or field components defined on a staggered 3D Yee grid.
    ftype : string
        If ftype=="scalar", values are considered to be scalars equally
        treated on the staggered grid. Interpolated value is the mean value of
        nearest next 12 neighbours -> from the center of each cell.
        If ftype=="vector", values are considered to be x,y,z components on the staggered grid.
        Interpolated vector contains components, where each one is the mean value of all next 4 neighbours
        of the corresponding component -> from the center of each cell.

    Returns
    -------
    f_intp : array
            Interpolated value/field data in 3D space
            
    """

    if arr.ndim != 4 or arr.shape[0] !=3:
        raise ValueError("Array must be 4-dimensional (number_field_comp, dim_x, dim_y, dim_z),\
                        where number_field_comp=3 and field values assumed on be distributed on \
                        a staggered Yee grid.")

    if not (ftype=="scalar" or ftype=="vector"):
        raise ValueError("ftype must be either 'scalar' or 'vector'.")

    fx, fy, fz = arr

    # Average each component over the next neighbours in each yee cell
    fx_interp = (fx[:-1, :-1, :-1] + fx[:-1, 1:, :-1] + fx[:-1, :-1, 1:] + fx[:-1, 1:, 1:]) / 4
    fy_interp = (fy[:-1, :-1, :-1] + fy[1:, :-1, :-1] + fy[:-1, :-1, 1:] + fy[1:, :-1, 1:]) / 4
    fz_interp = (fz[:-1, :-1, :-1] + fz[:-1, 1:, :-1] + fz[1:, :-1, :-1] + fz[1:, 1:, :-1]) / 4

    # Stack the arrays and compute the mean across the first axis
    f_comp_intp = np.stack([fx_interp, fy_interp, fz_interp])

    if (ftype == "scalar"):
        f_intp = np.expand_dims(np.mean(f_comp_intp, axis=0), 0)
    elif (ftype == "vector"):
        f_intp = f_comp_intp

    return f_intp


def calc_dissip_per_comp(wavelenght_nm_str, E, dens, params_ccpr):
    """
    Computation of the time-averaged electric power dissipation density 
    based on the electric field data from the FDTD simulation and the CCPR parameters.
    """  
    # Ensure that wavelenght is really given in units [nm] 
    matches = re.findall(r'(\d+)(\D+)', wavelenght_nm_str)
    wl, unit = matches[0] if matches else (None, None)
    if unit != "nm":
        raise ValueError("First argument expects wavelenght in nm followed \
                        by 'nm' as string, i.e. '<wavelenght>nm'")
    print("Computing dissipation per component for " + str(wl) + unit)
    wl = float(wl)
    # Calculate angular frequency (in SI units)
    omega = 2 * np.pi * V_LIGHT / (wl * 1e-9) 
    a = np.expand_dims(params_ccpr['a'], (1,2,3))
    c = np.expand_dims(params_ccpr['c'], (1,2,3))
    sigm = np.expand_dims(params_ccpr['sigm'], (1,2,3))
    eps = np.expand_dims(params_ccpr['eps'], (1,2,3))

    # Calculate dissipation, background material is assumed to be air (lossless) here
    dissip_per_comp = ((dens * 0.5 * sigm) + dens * np.sum(np.real((EPS_0 * c * omega**2) / ((1j*omega - a) * (-1j*omega - a))), axis=-1)) * np.abs(E)**2 
    
    return dissip_per_comp

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir_design', required=True, type=str, help='Directory from which the design data is read in')
    parser.add_argument('--store_q', required=False, action='store_true', help='Interpolate dissipation components and store the dissipation as .nc file')
    parser.add_argument('--store_dens', required=False, action='store_true', help='Interpolate density components and store the density as .nc file')
    parser.add_argument('--store_Eabs', required=False, action='store_true', help='Interpolate E field components and store the vector norm as .nc file')
    args = parser.parse_args()
    # ----------------------------------------------------
    # Read in CCPR parameters 
    p_ccpr = IO.readin_CCPR_params(os.path.join(args.dir_design, FILE_CCPR))
    
    # ----------------------------------------------------
    # Read in density and electric DFT fields
    dens = IO.readin_fields(args.dir_design, "dens")
    lambd, E_dict = IO.readin_fields(args.dir_design, "E") # lambd in m, E in V/m
    
    # ----------------------------------------------------
    # Reassign E fields, defined only in the design (& observation) region only
    params_sim = IO.readin_opt_regions(args.dir_design)
    drb, orb, or_in_drb = topt_regions_idx(params_sim) 
    
    # ----------------------------------------------------
    # Interpolate and store spatial density profile
    # ----------------------------------------------------
    if (args.store_dens):
        print("Interpolate density and store as NC file")
        dens_interpol_tmp = interpolate_components(dens, "scalar")
        IO.to_NetCDF(args.dir_design, dens_interpol_tmp, "dens_interpol.nc")   
    # ----------------------------------------------------
    # Crop E field array to the design region boundaries 
    # ---------------------------------------------------- 
    for w in E_dict:
        print("Cropping E " + w + " to DR boundaries")
        E_dict[w] = E_dict[w][:,drb[0,0]:drb[0,1] + 1, drb[1,0]:drb[1,1] + 1, drb[2,0]:drb[2,1] + 1]

        # Interpolate and store Eabs fields
        if (args.store_Eabs):
            print("Interpolate E vector" + w + " and store its norm as NC file")
            Eabs_interpol_tmp = np.linalg.norm(interpolate_components(E_dict[w], "vector"), axis=0, keepdims=True)
            IO.to_NetCDF(os.path.join(args.dir_design, "DFT"), Eabs_interpol_tmp, "Eabs_interpol_" + w + ".nc")
    
    # ----------------------------------------------------
    # Define integration parameters
    # ----------------------------------------------------
    Dx = params_sim['stepsizes'][0]
    Dy = params_sim['stepsizes'][1] 
    Dz = params_sim['stepsizes'][2]
    dxdydz = Dx * Dy * Dz
    # ----------------------------------------------------
    # Dissipation computation
    # ----------------------------------------------------
    diss_per_comp = {key: calc_dissip_per_comp(key, E_dict[key], dens, p_ccpr) for key in E_dict}
    if (args.store_q):
        diss_interpol = {key: np.expand_dims(np.sum(interpolate_components(diss_per_comp[key], "vector"), axis=0), 0) for key in diss_per_comp}
        for key, arr in diss_interpol.items():
            IO.to_NetCDF(os.path.join(args.dir_design, "DFT"), arr, "Dissip_interpol_" + key + ".nc")
    
    # Compute total dissipation by integration
    diss_integrated = {key: None for key in diss_per_comp}
    for key, arr in diss_per_comp.items():
        qx = np.sum(arr[0,:-1,:,:]) * dxdydz 
        qy = np.sum(arr[1,:,:,:])   * dxdydz  
        qz = np.sum(arr[2,:,:,:-1]) * dxdydz 
        q = qx + qy + qz
        diss_integrated[key] = q  
    # ----------------------------------------------------
    # Calculate, visualize and save Q_abs    
    # ----------------------------------------------------
    Wabs_integr = np.array(list(diss_integrated.values()))
    # get radius from file for normalization
    radius = IO.readin_radius(os.path.join(args.dir_design, FILE_SIM))
    prefac = np.pi * radius**2 * 0.5 * EPS_0 * V_LIGHT
    qabs_integr = Wabs_integr / prefac
    np.save(args.dir_design + "/lamb_qabs_radius", np.array([lambd, qabs_integr, radius], dtype=object))
    visualize.plt_quantity(args.dir_design, lambd, qabs_integr, "Qabs")

if __name__ == "__main__":
    main()