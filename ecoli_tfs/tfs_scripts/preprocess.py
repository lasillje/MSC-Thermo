"""Preprocessing utilities for preparing a COBRA model and thermodynamic inputs.

Contains:
  - Preprocess: class that loads model and external files (flux bounds,
    log-concentration priors, drG0 means/covariances), removes blocked
    reactions and unused metabolites, and produces constrained stoichiometric
    matrices and arrays used by downstream components.

The implementation mirrors the logic extracted from the original helper_pta.py
while being organized for testing and reuse.
"""

import os
import pandas as pd
import numpy as np
from cobra.flux_analysis.variability import find_blocked_reactions
from cobra.util.array import create_stoichiometric_matrix
from cobra.io import read_sbml_model
from pta.utils import ceil, floor

class Preprocess:
    """Loads files, sanitizes a COBRA model and provides prepared arrays.

    Parameters
    ----------
    cobra_file : str
        Path to SBML model file.
    vbound_file : str
        CSV with flux bounds indexed by reaction id, columns ['lb','ub'].
    logconcmean_file : str
        CSV with log-concentration mean estimates (EP output).
    logconcvar_file : str
        CSV with log-concentration variances (EP output).
    lncbounds_file : str
        CSV of TVA bounds for log-concentrations (used to detect fixed concs).
    drG0file : str
        CSV containing DrG'° means (equilibrator).
    drGcovsqrtfile : str
        CSV containing square-root of DrG'° covariance matrix.
    drG0covfile : str, optional
        Full covariance (not square-root) used optionally by EP-based parts.
    """

    def __init__(
        self,
        cobra_file,
        vbound_file,
        logconcmean_file,
        logconcvar_file,
        lncbounds_file,
        drG0file,
        drGcovsqrtfile,
        restrained_rxns=None,
        drG0covfile=None,
    ):
        if vbound_file is not None:
            data_bounds_flux = pd.read_csv(vbound_file, index_col=0)

        self.model = read_sbml_model(cobra_file)
        self.model.objective = self.model.reactions.biomass_EX

        for x in self.model.reactions:
            print(x.id, x.lower_bound, x.upper_bound)

        # Load drG0 info from provided files
        drg0_prime_mean_init = pd.read_csv(drG0file, index_col=0).values.T[0]
        drg0_prime_cov_sqrt_init = pd.read_csv(
            drGcovsqrtfile, index_col=0
        ).values
        if drG0covfile is not None:
            drg0_prime_cov_init = pd.read_csv(drG0covfile, index_col=0).values

        print(len(drg0_prime_mean_init))
        # Record initial metabolites and remove specific ions (Mg)
        initial_mets_id = [m.id for m in self.model.metabolites]
        # remove Mg if present (safe)
        for mg_id in ("Mg_c", "Mg_e", "Mg_m"):
            if mg_id in [m.id for m in self.model.metabolites]:
                try:
                    self.model.metabolites.remove(self.model.metabolites.get_by_id(mg_id))
                except Exception:
                    pass

        # Adjust model bounds to provided flux estimates and add margins
        #self.add_v_margin(data_bounds_flux)

        #self.model.reactions.get_by_id("ADPATPtm-4-2").bounds = (0.0, 0.0)
        #print('Adjusting bounds for ADPATPtm-4-2',self.model.reactions.get_by_id("ADPATPtm-4-2").bounds)
        #self.model.reactions.get_by_id('G3PT').bounds=(0,20)
        #self.model.reactions.get_by_id('GCY1').bounds=(0,20)
        
        # Find blocked reactions and remove them (remember original indices)
        init_reactions_id = [r.id for r in self.model.reactions]
        id_removed = find_blocked_reactions(self.model)
        print(id_removed)
        rxn_removed = {
            old_idx: rid for old_idx, rid in enumerate(init_reactions_id) if rid in id_removed
        }
        ridx_notremoved = [
            oldidx for oldidx in range(len(init_reactions_id)) if oldidx not in rxn_removed.keys()
        ]
        old_new_r_mapping = {old: new for new, old in enumerate(ridx_notremoved)}
        # Remove blocked reactions and orphan metabolites
        self.model.remove_reactions(id_removed, remove_orphans=True)

        # Track metabolites removed compared to the initial list
        mets_removed = {
            old_idx: mid
            for old_idx, mid in enumerate(initial_mets_id)
            if mid not in [m.id for m in self.model.metabolites]
        }
        met_idx_notremoved = [
            oldidx for oldidx in range(len(initial_mets_id)) if oldidx not in mets_removed.keys()
        ]

        # Stoichiometric matrix for the cleaned model
        S = create_stoichiometric_matrix(self.model)
        self.S = S

        # Reactions for which the second law is ignored (original indices)
        # Restrict to FCA reactions
        # ~45 points
        if restrained_rxns is not None:
            snd_ignored_idxs = restrained_rxns
        else:
            snd_ignored_idxs = ['biomass_EX', 'EX_o2', 'biomass_ce', 'EX_ac', 'biomass', 'EX_so4', 'EX_oro', 'H2Ot', 'EX_pi', 'EX_nh3', 'EX_co2', 'EX_h', 'EX_glc', 'EX_h2o']
        
        #snd_not_involved = [x.id for x in self.model.reactions if x.upper_bound > 90 or x.lower_bound < -90]
        
        #print(snd_not_involved)

        #snd_ignored_idxs += snd_not_involved
        
        print(snd_ignored_idxs)
        print(len(snd_ignored_idxs))
        
        #snd_ignored_idxs = [x.id for x in self.model.reactions if "biomass" in x.id or "EX_" in x.id or "BIOMASS" in x.id]
        #snd_ignored_idxs = ['EX_h2o_e', 'EX_nh4_e', 'EX_glu__L_e', 'EX_lac__D_e', 'EX_acald_e', 'EX_akg_e', 'biomass_ce', 'EX_for_e', 'EX_ac_e', 'EX_glc__D_e', 'BIOMASS_Ecoli_core_w_GAM', 'EX_pyr_e', 'EX_succ_e', 'biomass_EX', 'EX_etoh_e', 'H2Ot', 'EX_pi_e', 'EX_co2_e', 'EX_h_e', 'EX_o2_e']

        snd_ignored = [ self.model.reactions.index(self.model.reactions.get_by_id(x)) for x in snd_ignored_idxs ]
        print(snd_ignored)
        #pos_ht_rxn = [205, 206, 207, 208] ##posqt and HT are fully determined by drG0prime 
        all_ignored = snd_ignored #+ pos_ht_rxn
        snd_ignored_newidx = [
            old_new_r_mapping[old] for old in all_ignored if old in ridx_notremoved
        ]

        # Constrained reaction indices and mapping
        rxn_constrained_idx = [ridx for ridx in range(len(self.model.reactions)) if ridx not in snd_ignored_newidx]
        self.Sconstrained = S[:, [idx for idx in range(len(self.model.reactions)) if idx not in snd_ignored_newidx]]
        self.rid_constrained = [self.model.reactions[idx].id for idx in range(len(self.model.reactions)) if idx not in snd_ignored_newidx]

        # ---- Load concentration priors ----

        if lncbounds_file is not None:
            bounds_lnconc = pd.read_csv(lncbounds_file, index_col=0).loc[[m.id for m in self.model.metabolites]]
            ignore_met_id = bounds_lnconc[bounds_lnconc["lb"] > 900].index
            ignoreconc_met_idx = [self.model.metabolites.index(met_id) for met_id in ignore_met_id]
            self.ignoreconc_met_idx = ignoreconc_met_idx


        log_conc_mean = pd.read_csv(logconcmean_file, index_col=0)["0"].values
        log_conc_var = pd.read_csv(logconcvar_file, index_col=0)["0"].values

        log_conc_mean = [x for i, x in enumerate(log_conc_mean) if i not in self.ignoreconc_met_idx]
        log_conc_var = [x for i, x in enumerate(log_conc_var) if i not in self.ignoreconc_met_idx]

        log_conc_cov = np.diag(log_conc_var)
        
        self.log_conc_cov = log_conc_cov
        self.log_conc_mean = log_conc_mean

        print(len(self.log_conc_mean), self.log_conc_mean)
        print(len(self.log_conc_cov), self.log_conc_cov)

        # ---- Load drG'° statistics restricted to constrained reactions ----
        idx_constrained_old = [
            idx for idx in range(len(drg0_prime_mean_init)) if idx not in all_ignored and idx in ridx_notremoved
        ]

        print(idx_constrained_old)

        self.drg0_prime_mean = drg0_prime_mean_init[idx_constrained_old]
        self.drg0_prime_cov_sqrt = drg0_prime_cov_sqrt_init[idx_constrained_old, :]#[:, idx_constrained_old]
        if drG0covfile is not None:
            self.drG0_prime_cov = drg0_prime_cov_init[idx_constrained_old, :][:, idx_constrained_old]

    def add_v_margin(self, data_bounds_flux, margin: float = 1e-4, round_to_digits: int = 6):
        """Apply small margins to bounds read from a CSV to avoid exact-tight bounds.

        This mirrors the behaviour in the original monolithic script to avoid
        numerical edge cases.
        """
        for r in self.model.reactions:
            lb = data_bounds_flux.loc[r.id, "lb"]
            ub = data_bounds_flux.loc[r.id, "ub"]

            if 0 <= lb <= margin:
                r.lower_bound = max(r.lower_bound, 0.0)
            else:
                r.lower_bound = max(r.lower_bound, floor(lb - margin, round_to_digits))

            if -margin <= ub <= 0:
                r.upper_bound = min(r.upper_bound, 0.0)
            else:
                r.upper_bound = min(r.upper_bound, ceil(ub + margin, round_to_digits))

