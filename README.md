# Documentation for dissipation-calculation {#README}

This code allows for calculating and visualizing the electric power dissipation profile and its integration (absorption efficiency) of the TopOpt designs and spheres made of Gold and Silicon presented in the article Ref.XXX. The computations require the dataset provided at XXX, containing the electric field distribution in 3D for multiple wavelengths (netCDF), the final density (netCDF), the design (STL), and the material and simulation parameters (JSON) used in the optimization.

## Requirements 
- Download this repository as well as the dataset (under "Data and Resources") from the following URL: XXX
- Replace the "data" directory with the downloaded directory sharing the same name.
- Create an environment (e.g., conda environment) with the Python version and necessary packages that are listed in `requirements.txt`.

## Usage
The file `commands.txt` provides examples of prepared commands one can execute:
- `analyze.py` computes and visualizes the absorption efficiency $Q_{\text{abs}}$. In addition, the flags `--store_q`, `--store_Eabs`, and `--store_dens` can be used to store the spatially interpolated data (based on the original one defined on the staggered Yee grid) for further computations or visualization as NetCDF files.
- `verifyEpsilon.py` serves as a verification of the permittivity fit using the complex-conjugate pole-residue (CCPR) model.
- `verifyMie.py` is used exclusively to verify the comparison of the absorption efficiency of a sphere with Mie theory. Make sure that you have first run `analyze.py` for the data set of the sphere, as the stored output is required for this calculation specifically.

## Examples
### TopOpt design (Gold)
Run: ``python analyze.py --dir_design=data/Au/Design``
<p align="left">
<img src="https://github.com/JoGed/dissipation-calculation/assets/83292544/b9971998-b672-4f53-90c7-9ab63489d96d" alt="Qabs_design_gold" width="300">
</p>
<p>Output: Absorption efficiency $Q_{\text{abs}}$ of the design.</p>

<br>

Run: ``python verifyEpsilon.py --lower_bound=0.3 --upper_bound=0.7 --dir_CCPR=data/Au/Design/materialTopOpt.json --dir_data=data/eps_data/Au_Johnson.txt --N=201``
<p align="left">
<img src="https://github.com/JoGed/dissipation-calculation/assets/83292544/d4929d3c-2bac-40f0-aabd-739427f56100" alt="eps_gold" width="300">
</p>
<p>Output: Permittivity $\varepsilon$ fit of Gold.</p>

### Sphere (Gold)
Run: ``python analyze.py --dir_design=data/Au/Sphere``
<p align="left">
<img src="https://github.com/JoGed/dissipation-calculation/assets/83292544/4ff32c4f-461a-4637-a9df-64dd37ef9984" alt="Qabs_sphere_gold" width="300">
</p>
<p>Output: Absorption efficiency $Q_{\text{abs}}$ of a sphere ($d$=300nm).</p>

<br>

Run: ``python verifyMie.py --dir_mie=data/Au/Sphere``
<p align="left">
<img src="https://github.com/JoGed/dissipation-calculation/assets/83292544/316d7692-eb93-4031-9a66-c3c5c9c09134" alt="mie_gold" width="300">
</p>
<p>Output: Absorption efficiencies $Q_{\text{abs}}$ of the sphere compared to Mie theory ($d$=300nm).</p>

<br>

