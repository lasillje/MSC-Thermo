# MSC-Thermo

## Workflow

# Loading model

The base model is loaded using Thermo-Flux from the data in /datafiles/model.xlsx
Additional experimental and metabolome data is loaded from the files in /regression/
A comprehensive function to load a model based on various input files can be found in scripts/gen_model.py -> gen_model()

# Applying physiological data

Additional physiological data and metabolome data can be applied to the base model. This can be done using the function apply_physio_data() found in scripts/gen_model.py
This requires a functional thermo-flux model to already be loaded on which the physio and metabolome data will be applied

# Removing blocked reactions

After physiological and metabolome data is applied, blocked reactions are listed using list_blocked_reactions() from scripts/reaction_utils.py
This returns a list of blocked reactions, in this case we don't allow other excretions in our model, and we do not open exchange fluxes when finding blocked reactions.
After this, all blocked reactions are removed from the model, using tmodel.remove_reactions(remove_orphans=True). Note that remove_orphans must be True, otherwise the model becomes infeasible.
This is done to make the model as small but functional as possible, saving disk space and optimization time down the line. Another approach, setting the flux bounds of blocked reactions to [0, 0], was discarded as the optimization time would exceed 6+ hours for some reactions.