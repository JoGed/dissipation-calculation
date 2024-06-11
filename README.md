# Documentation for dissipation-calculation {#README}

This code allows to calculate and vizualise the electric power dissipation profile and its integration (absorption effiency)
of the TopOpt designs + spheres made of Gold and Silicon presented in the article Ref.XXX.
The compuations require the dataset provided at XXX, containing the electric field distribution in 3D for multiple wavelengths (netCDF), 
the final density (netCDF), the design (STL) and material and simulation parameters (JSON) used in the optimization.

## Requirements {#requirements}
- Download this repository as well the dataset (under "Data and Resources") from the following URL:
- Replace the "data" dircetory with the downloaded directory sharing the same name
- Create an environment (e.g. conda environment) with the python version and necessary packages that are listed in `requirements.txt`

## Usage
The file `commands.txt` provide examples of preprepared commands one can execute
- `analyze.py`, computes and vizualises the absorption effiency $Q_{\text{abs}}$. In addition, the flags `--store_q`, `--store_Eabs`, and `--store_dens`can be used to store the spatially interpolated data (based on the original one defined on the staggered Yee grid) for further computations or vizualisation.
- `verifyEpsilon.py` serves as a verification of the permittivity fit using the complex-conjugate pole-residue (CCPR) model
- `verifyMie.py` is used exclusively to verify the comparison of the absoprtion efficiency of a sphere with Mie theory. Make sure that you have first run analyze.py for the data set of a simple sphere, as the numerical data is required for the calculation.
### Example
