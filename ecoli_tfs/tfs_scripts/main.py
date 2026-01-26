"""Orchestration: build a SupTFSmodel from input files.

Example
-------
from helper_pta import main
tfs_model = main(cobrafile, vbounds_file, logconcmean_file, logconcvar_file,
                 lncbounds_file, drG0file, drGcovfile, ranktol=1e-5)

The function returns an initialized SupTFSmodel instance (not yet sampled).
"""

import pta
import cvxpy as cp
from .preprocess import Preprocess
from .thermospace_mod import ThermodynamicSpaceMod
from .tbasis_mod import ThermodynamicSpaceBasismod
from .tfsmodel_mod import SupTFSmodel
from .pmo_mod import PmoProblemMod
#from sampler import helper_sampling as hs

def main(cobrafile, vbounds_file,
         logconcmean_file, logconcvar_file, lncbounds_file,
         drG0file, drGcovfile, ranktol=1e-5, restrained_rxns=None, sector_dict=None):
    """Create and return a configured SupTFSmodel.

    Parameters match those of the original script. The returned object is ready
    for sampling (call .sample(points_array) on it), but no sampling is run by
    this function.
    """
    prep = Preprocess(
        cobrafile,
        vbounds_file,
        logconcmean_file,
        logconcvar_file,
        lncbounds_file,
        drG0file,
        drGcovfile,
        restrained_rxns=restrained_rxns
    )
    # Check TFVA bounds applied
    # Ignore conc/met important

    #if sector_dict is not None:
    #   hs.set_sector_flux(prep.model, sector_dict)

    thermodynamic_space = ThermodynamicSpaceMod(
        prep.Sconstrained,
        prep.rid_constrained,
        prep.model.metabolites,
        prep.drg0_prime_mean,
        prep.drg0_prime_cov_sqrt,
        prep.log_conc_cov,
        prep.log_conc_mean,
    )

    tbasis = ThermodynamicSpaceBasismod(
        thermodynamic_space,
        explicit_log_conc=False,
        explicit_drg0=False,
        explicit_drg=True,
        min_eigenvalue=1e-10,
        ranktol=ranktol,
        ignoreconc_met_idx=prep.ignoreconc_met_idx,
    )

    problem = PmoProblemMod(prep.model, thermodynamic_space, tbasis, solver="GUROBI")

    tfs_model = SupTFSmodel(prep.model, thermodynamic_space, tbasis, problem, solver=cp.GUROBI)
    tfs_model.ignoreconc_met_idx = prep.ignoreconc_met_idx

    return tfs_model
