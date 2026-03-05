# MSC-Thermo

## Workflow

The repository consists of two main jupyter notebooks, present in the ecoli_flux_bounds and ecoli_tfs directories. The intended progression is to first go through the ecoli_flux_bounds notebook to subsequently use the generated data in the ecoli_tfs notebook.

### Thermodynamically constrained variability analyses
The notebook within the ecoli_flux_bounds directory first provides a walkthrough of setting up a stoichiometric-thermodynamic model using Thermo-Flux, using the various datafiles found in ecoli_flux_bounds/datafiles/. Furthermore, this notebook contains code to perform thermodynamically constrained variability analyses (TVA) of metabolite concentrations and reaction fluxes (TFVA). The first part of the notebook sets up a stoichiometric-thermodynamic E. coli model using metabolome data and experimental flux data for two different physiological growth conditions (glucose and acetate). After this, for both conditions, metabolite TVA and TFVA is setup and performed using Thermo-Flux functionality. Lastly, a graphing section is present at the end that combines the results from the metabolite TVA and TFVA into a figure. The ecoli_flux_bounds/results/ folder holds the raw output data from the metabolite TVA and TFVA for both growth conditions (*_objval.txt files). Furthermore, the ecoli_flux_bounds/results/bounds/ directory holds the intermediate further processed results, with the ecoli_flux_bounds/results/bounds/composite/ directory holding the final TFVA data files for both conditions, currently these are the *_final_tfva.csv files, which are made through the combination of the *_metabolome_tfva.csv and *_no_metabolome_tfva.csv files (saved within the notebook using the merge_width_no_metabolome function).


### Thermodynamically constrained flux sampling
The notebook within the ecoli_tfs directory builds upon the previous notebook, using the metabolite TVA and TFVA results to perform thermodynamically constrained flux sampling using the PTA library (installation instructions below). Stoichiometric-thermodynamic models are again setup using the ecoli_tfs/datafiles/ files, in addition to using the completed metabolite TVA and TFVA files in ecoli_tfs/results and ecoli_tfs/results/bounds respectively. Thermodynamically constrained flux sampling is setup and performed using the PTA library, of which the notebook provides a walkthrough of generating the necessary files to perform this analysis in a timely manner, as the used model is quite large. To this end, additional flux-coupling analysis is also performed, which provides a network of fully coupled reactions (flux of reaction A has constant effect on flux of reaction B), this way, the total amount of reactiosn that need to be thermodynamically constrained in the flux sampling can be reduced, speeding up the process. Lastly, a graphing part is added as well using the resulting thermodynamically constrained flux sampling distributions. Files that I used on the cluster are also provided (tfs_points.py, tfs_drg_and_fluxes.py), as well as their respective slurm script to set up the jobs on the cluster.

# PTA installation

For thermodynamically constrained flux sampling, the PTA package is used (https://gitlab.com/csb.ethz/pta). Installation on Windows is not possible, so on Windows I used WSL. Then, manually build the PTA package from the source.
With this, issues may still arise. The exact steps I followed to get it to work:
First, if not on the cluster then install the gurobi python interface with 
```bash
conda install -c gurobi gurobi
```

Otherwise, on the cluster you can use
```bash
ml Gurobi
ml Ninja
```

Then, download the gurobi c/c++ files from the gurobi website (linux version, .tar.gz) and unzip it

Then, clone the PTA repository and set the CMAKE minimum version in pta/cpp/python/CMakeLists.txt to 3.5

Make sure to also update and load all submodules while in the pta/ directory
```bash
git submodule update --init --recursive
```

Next, make sure to run
 ```bash
conda install wheel setuptools
``` 
if not yet installed. 

Next, In the pta/cpp/python/pybind11 we need to update the pybind version:
```bash
git fetch --all
git checkout v2.11.1
```
Then we need to tell PTA where our gurobi source download is, for example:

```bash
export GUROBI_HOME=/home/user/gurobi/gurobi1300/linux64
```

This should also be put in either ~/.profile or ~/.bashrc

After this, building needs to be ran in isolated mode and the build should succeed:
```bash
CMAKE_PREFIX_PATH=$GUROBI_HOME pip install . --no-build-isolation
```

Now stuff like
```python
import enkie
import pta

model = pta.load_example_model("e_coli_core")
model.reactions.BIOMASS_Ecoli_core_w_GAM.lower_bound = 0.5
```

should work.

However, some functions also require Java to be installed, so if not done yet:

```bash
conda install openjdk
java -version
```

and add the following to ~/.profile or ~/.bashrc :
```bash
export JAVA_HOME=$CONDA_PREFIX/lib/jvm
```

After all this hassle PTA should work.
