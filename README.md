# MSC-Thermo

## Workflow

# Loading model

The base model is loaded using Thermo-Flux from the data in /datafiles/model.xlsx
Additional experimental and metabolome data is loaded from the files in /regression/
A comprehensive function to load a model based on various input files can be found in scripts/gen_model.py -> gen_model()

# Applying physiological data

Additional physiological data and metabolome data can be applied to the base model. This can be done using the function apply_physio_data() found in scripts/gen_model.py
This requires a functional thermo-flux model to already be loaded on which the physio and metabolome data will be applied

This is done for various growth conditions: WT-Glc_I, WT-Ace_I

# Removing blocked reactions

After physiological and metabolome data is applied, blocked reactions are listed using list_blocked_reactions() from scripts/reaction_utils.py
This returns a list of blocked reactions, in this case we don't allow other excretions in our model, and we do not open exchange fluxes when finding blocked reactions.
After this, all blocked reactions are removed from the model, using tmodel.remove_reactions(remove_orphans=True). Note that remove_orphans must be True, otherwise the model becomes infeasible.
This is done to make the model as small but functional as possible, saving disk space and optimization time down the line. Another approach, setting the flux bounds of blocked reactions to [0, 0], was discarded as the optimization time would exceed 6+ hours for some reactions.

# Thermodynamic FVA

For each growth condition, thermodynamic FVA (TFVA) was done to constrain flux space further. For each condition, the individual reactions are output as gurobi model with two scenarios (minimize and maximize). This is done in analysis_thermo_fva.ipynb. The resulting gurobi model files are uploaded to the cluster on which they are optimized. After this, bounds are updated for each reaction in the thermo-flux model using tfva_update_bounds() from scripts/reaction_utils

# Thermodynamic Sampling

For thermodynamic smapling, the PTA package is used (https://gitlab.com/csb.ethz/pta). Installation on Windows is not possible, so on Windows I used WSL. Then, manually build the PTA package from the source.
With this, issues may still arise. The exact steps I followed to get it to work:
First, install the gurobi python interface with 
```bash
conda install -c gurobi gurobi
```
Then, download the gurobi c/c++ files from the gurobi website (linux version, .tar.gz) and unzip it
Then, clone the PTA repository and set the CMAKE minimum version in pta/cpp/python/CMakeLists.txt to 3.5
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
