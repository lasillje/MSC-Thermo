"""Thermodynamic-space basis computation (modified).

This module contains `ThermodynamicSpaceBasismod`, a modified version of the
PTA basis computation that supports ignoring a set of metabolite indices
(useful when some concentrations are fixed by TVA).
"""

import numpy as np
import scipy.linalg
import numpy.linalg as la
from pta.utils import covariance_square_root
from component_contribution.linalg import LINALG
from enkie.constants import R
import pta

class ThermodynamicSpaceBasismod(pta.ThermodynamicSpaceBasis):
    """Compute a minimal QR basis for the (log_conc, drg0, drg) observables.

    Parameters
    ----------
    thermodynamic_space : ThermodynamicSpaceMod
        Input container built by ThermodynamicSpaceMod.
    explicit_log_conc : bool
    explicit_drg0 : bool
    explicit_drg : bool
    min_eigenvalue : float
        Tolerance passed to covariance_square_root.
    ranktol : float
        Tolerance used by the LINALG.qr_rank_deficient routine.
    ignoreconc_met_idx : list
        Indices of metabolites to ignore when building S (their columns are removed).
    """

    def __init__(
        self,
        thermodynamic_space,
        explicit_log_conc=True,
        explicit_drg0=True,
        explicit_drg=True,
        min_eigenvalue: float = 1e-10,
        ranktol: float = 1e-5,
        ignoreconc_met_idx: list = None,
    ):
        if ignoreconc_met_idx is None:
            ignoreconc_met_idx = []

        self._to_log_conc_transform = None
        self._to_drg0_transform = None
        self._to_drg_transform = None

        self._observables_ranges = {"drg": [], "drg0": [], "log_conc": []}
        self.ignoreconc_met_idx = ignoreconc_met_idx

        # compute the basis (mirrors the logic from the monolithic file)
        self._compute_basis_mod(
            thermodynamic_space,
            explicit_log_conc,
            explicit_drg0,
            explicit_drg,
            min_eigenvalue,
            ranktol,
        )

    def _compute_basis_mod(
        self,
        thermodynamic_space,
        explicit_log_conc,
        explicit_drg0,
        explicit_drg,
        tolerance,
        ranktol,
        T=303.15,
    ):
        # Extract numerical arrays
        log_conc_mean = thermodynamic_space.log_conc_mean
        log_conc_cov = thermodynamic_space.log_conc_cov
        drg0_prime_mean = thermodynamic_space.drg0_prime_mean
        drg0_prime_cov_sqrt = thermodynamic_space.drg0_prime_cov_sqrt
        S_con_T = thermodynamic_space.S_constraints.T

        # If there are ignored metabolite indices, remove their columns from S
        if len(self.ignoreconc_met_idx) > 0:
            keep_idxs = [idx for idx in range(S_con_T.shape[1]) if idx not in self.ignoreconc_met_idx]
            S_con_T_mod = S_con_T[:, keep_idxs]
        else:
            S_con_T_mod = S_con_T

        RT = (R.m * T)
        n_mets = log_conc_mean.size
        n_rxns = drg0_prime_mean.size

        # Compute square-root covariance for log concentrations
        log_conc_cov_sqrt = covariance_square_root(log_conc_cov, tolerance)

        # Combine means and square-root covariance
        log_conc_drg0_mean = np.block([[log_conc_mean], [drg0_prime_mean]])
        log_conc_drg0_cov_sqrt = scipy.linalg.block_diag(log_conc_cov_sqrt, drg0_prime_cov_sqrt)

        print(log_conc_drg0_mean)
    
        # Build mapping to observables
        mapping_blocks = []
        if explicit_log_conc:
            mapping_blocks.append([np.identity(n_mets), np.zeros((n_mets, n_rxns))])
        if explicit_drg0:
            mapping_blocks.append([np.zeros((n_rxns, n_mets)), np.identity(n_rxns)])
        if explicit_drg:
            mapping_blocks.append([RT * S_con_T_mod, np.identity(n_rxns)])

        to_observables = np.block(mapping_blocks)
        observables_mean = to_observables @ log_conc_drg0_mean
        observables_cov_sqrt = to_observables @ log_conc_drg0_cov_sqrt

        # Minimal square-root via rank-deficient QR from component_contribution
        observables_cov_sqrt = LINALG.qr_rank_deficient(observables_cov_sqrt.T, ranktol).T
        observables_cov_sqrt[np.abs(observables_cov_sqrt) < 1e-5] = 0

        self._to_observables_transform = (observables_cov_sqrt, observables_mean)
        self._dimensionality = observables_cov_sqrt.shape[1]
        self._sigmas = la.norm(observables_cov_sqrt, axis=0).T

        # expose helpful attributes for downstream EP/PMO code
        self.to_observables = to_observables
        self.log_conc_drg0_mean = log_conc_drg0_mean
        self.log_conc_drg0_cov_sqrt = log_conc_drg0_cov_sqrt

        # Build transforms to each observable block
        current_id = 0
        if explicit_log_conc:
            self._observables_ranges["log_conc"] = list(range(current_id, current_id + n_mets))
            self._to_log_conc_transform = (
                observables_cov_sqrt[self._observables_ranges["log_conc"], :],
                observables_mean[self._observables_ranges["log_conc"]],
            )
            current_id += n_mets
        if explicit_drg0:
            self._observables_ranges["drg0"] = list(range(current_id, current_id + n_rxns))
            self._to_drg0_transform = (
                observables_cov_sqrt[self._observables_ranges["drg0"], :],
                observables_mean[self._observables_ranges["drg0"]],
            )
            current_id += n_rxns
        if explicit_drg:
            self._observables_ranges["drg"] = list(range(current_id, current_id + n_rxns))
            self._to_drg_transform = (
                observables_cov_sqrt[self._observables_ranges["drg"], :],
                observables_mean[self._observables_ranges["drg"]],
            )
