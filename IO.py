import netCDF4 as nc
import numpy as np
import json
import glob
import os
import re
from constants import FILE_DENS, FILE_DFTS, FILE_SIM

def readin_CCPR_params(directory):
    """
      Returns dictionary with the CCPR parameters read from a json file
    """
    file_str = directory 
    # Open the file and read the JSON data
    with open(file_str, 'r') as file:
      material_dict = json.load(file)[0]
    s = ["x", "y", "z"]
    eps = np.array([material_dict[f'eps{axis}1'] for axis in s])
    sigm = np.array([material_dict[f'sigm{axis}1'] for axis in s])
    damp = np.array([material_dict[f'damp{axis}1'] for axis in s])
    cr = np.array([material_dict[f'c{axis}r1'] for axis in s])
    ci = np.array([material_dict[f'c{axis}i1'] for axis in s])
    ar = np.array([material_dict[f'a{axis}r1'] for axis in s])
    ai = np.array([material_dict[f'a{axis}i1'] for axis in s])
    c = cr + 1j*ci
    a = ar + 1j*ai

    dict = {
        'eps': eps,
        'sigm': sigm,
        'damp': damp,
        'c': c,
        'a': a
    }

    return dict


def readin_fields(directory, field_type):
  """
    Return wavelength and fields (as dictionary) based on the 3D electric fields
    in the frequency domain, read in from the dataset stored in the directory
  """
  if (field_type == "dens" or field_type == "gradient"):
    # Find files that match the pattern
    print("Readin densitiy")
    dens_file_path  = glob.glob(os.path.join(directory, FILE_DENS))[0]
    dataset = nc.Dataset(dens_file_path)
    fields = np.array([dataset.variables['Ex'], 
                     dataset.variables['Ey'], 
                     dataset.variables['Ez']
            ])   
    return fields
    
  elif (field_type == "E"):

    # Find files that match the pattern 
    E_file_paths = glob.glob(os.path.join(directory, FILE_DFTS))
    wavelength_str = [os.path.basename(file_path).split('_')[1] for file_path in E_file_paths]
    E = []
    wavelength = []

    for i in range(len(E_file_paths)):
      matches = re.findall(r'(\d+)(\D+)', wavelength_str[i])
      wl, unit = matches[0] if matches else (None, None)
      if unit != "nm":
        raise ValueError("E field file must have specified name \
                          containing the wavelenght in nm, i.e. 'E_<wavelength>nm'")
      print("Readin E for " + str(wl) + unit)
      wavelength.append(float(wl))

      with nc.Dataset(E_file_paths[i]) as dataset:
          E_real = []
          E_imag = []
          for component in ['Ex', 'Ey', 'Ez']:
              E_real.append(dataset.variables[f'{component}_real'])
              E_imag.append(dataset.variables[f'{component}_imag'])

          E.append(np.array(E_real) + 1j * np.array(E_imag))

    wl_sorted, wl_str_sorted, E_sorted = zip(*sorted(zip(wavelength, wavelength_str, E)))
    
    print("Create sorted E dictionary...")
    fields = {}
    for i in range(len(wl_str_sorted)):
        key = wl_str_sorted[i]
        value = E_sorted[i]
        fields[key] = value

    # Convert to m instead of nm
    wl_sorted = np.array(wl_sorted) * 1e-9
    return wl_sorted, fields


def readin_opt_regions(directory):
  """
    Read in boundaries of design and observation region, 
    and domain (including CPML) within the FDTD setup
  """
  file_str = os.path.join(directory, FILE_SIM)
  with open(file_str, 'r') as file:
    params_dict = json.load(file)[0]

  s = ["X", "Y", "Z"]
  design_region = np.array([params_dict[f'DR 0 in {axis} [Idx]'] for axis in s])
  observation_region = np.array([params_dict[f'OR 0 in {axis} [Idx]'] for axis in s])
  fdtd_domain_size = np.array(params_dict['Domain Size'])
  cpml_size = np.array(params_dict['PML Box'])
  spatial_steps = np.array(params_dict['Space Step'])

  dict = {
        'dr': design_region,
        'or': observation_region,
        'domain': fdtd_domain_size,
        'cpml': cpml_size,
        'stepsizes': spatial_steps
        }
  
  return dict

def to_NetCDF(directory, arr, name_str):
  """
    Converts array to NetCDF file and stores it
  """
  if arr.ndim != 4:
    raise ValueError("Array must be 4-dimensional (number_field_comp, dim_x, dim_y, dim_z)")
   
  with nc.Dataset(os.path.join(directory, name_str), 'w', format='NETCDF4') as dataset:
    dim_var = ('x','y','z')
    # Create dimensions
    for dim_name, size in zip(dim_var, arr.shape[1:]):
      dataset.createDimension(dim_name, size)
    print("NetCDF output " + name_str)
    is_complex = np.issubdtype(arr.dtype, np.complexfloating)
    n_components = arr.shape[0]
    if (is_complex): float_type = arr[0, :, :, :].real.dtype

    for i in range(n_components):
        component_suffix = f'_{dim_var[i]}' if n_components > 1 else ''
        var_base_name = f'f{component_suffix}'

        if is_complex:
            # Create separate variables for real and imaginary parts
            dataset.createVariable(f'{var_base_name}_real', float_type, dim_var)[:, :, :] = np.real(arr[i, :, :, :])
            dataset.createVariable(f'{var_base_name}_imag', float_type, dim_var)[:, :, :] = np.imag(arr[i, :, :, :])
        else:
            # Create a single variable for real data
            dataset.createVariable(var_base_name, arr.dtype, dim_var)[:, :, :] = arr[i, :, :, :]
  return 1

def readin_radius(file_path):
  """
    Reads in the radius from a json file
  """
  with open(file_path, 'r') as file:
    params_dict = json.load(file)[0]
  radius = params_dict['Radius']
  return radius

def readin_wavelength(file_path):
  """
    Reads in the number of wavelengths and boundary values from a json file
  """
  with open(file_path, 'r') as file:
    params_dict = json.load(file)[0]
  wl_min = params_dict['Minimum Wavelength']	
  wl_max = params_dict['Maximum Wavelength']
  wl_N = 	int(params_dict['Number of Wavelengths'])	
  lambd = np.linspace(wl_min, wl_max, num=2*wl_N - 1 , endpoint=True)
  return lambd
   