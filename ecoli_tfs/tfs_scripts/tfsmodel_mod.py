"""Sampling wrapper (SupTFSmodel) that extends the PTA TFSModel.

This module provides SupTFSmodel which:
  - stores the original COBRA model reference,
  - prepares sampling settings,
  - finds reversible/irreversible direction candidates and
  - wraps the convergence manager logic.
"""

import numpy as np
import logging
from pta.sampling import tfs
import cvxpy as cp
from . import tbasis_mod
from . import conv_manager_mod as cm
from pta.constants import (
    default_max_psrf,
    default_max_threads,
    default_min_eigenvalue_tds_basis,
    tfs_default_feasibility_cache_size,
    tfs_default_min_rel_region_length,
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

class SupTFSmodel(tfs.TFSModel):
    """Extended TFSModel with conveniences for the clustered sampler.

    Parameters
    ----------
    model : cobra.Model
    thermodynamic_space : ThermodynamicSpaceMod
    tbasis : ThermodynamicSpaceBasismod
    problem : pta.PmoProblem
    Additional sampler parameters can be provided as keyword arguments.
    """
    
    def __init__(self, model, thermodynamic_space, tbasis, problem, num_samples=int(5e5),
                 num_dir_samples=int(5e7), num_initial_steps=None, max_steps=2**25,
                 max_threads=32, warmup=2, **kwargs):
        super().__init__(model, thermodynamic_space, tbasis, **kwargs)
        self.model = model
        self.problem = problem
        self.num_samples = num_samples
        self.num_dir_samples = num_dir_samples
        self.max_steps = max_steps
        self.num_initial_steps = num_initial_steps
        self.max_threads = max_threads
        self.warmup = warmup

        # prepare directions
        self.find_rev_irr_dirs()
        self.directions_transform = None

    def find_rev_irr_dirs(self):
        reaction_idxs_T = list(range(len(self.T.reaction_ids)))
        reaction_idxs_F = [self.F.reaction_ids.index(id) for id in self.T.reaction_ids]
        only_forward_ids_T = [i for i in reaction_idxs_T if self.F.lb[reaction_idxs_F[i]] >= 0]
        only_backward_ids_T = [i for i in reaction_idxs_T if self.F.ub[reaction_idxs_F[i]] <= 0]
        reversible_ids_T = [i for i in reaction_idxs_T if self.F.lb[reaction_idxs_F[i]] < 0 and self.F.ub[reaction_idxs_F[i]] > 0]

        reversible_dirs = [(i, -1) for i in reversible_ids_T] + [(i, 1) for i in reversible_ids_T]
        irreversible_dirs = [(i, -1) for i in only_backward_ids_T] + [(i, 1) for i in only_forward_ids_T]

        self.reversible_dirs = reversible_dirs
        self.irreversible_dirs = irreversible_dirs

    def safe_find_point(self, direction):
        try:
            return tfs._find_point(self.problem, direction)
        except Exception as e:
            logger.warning("safe_find_point failed for %s: %s", direction, e)
            return None

    def sampling_settings(self, points_array):
        settings = tfs.pb.FreeEnergySamplerSettings()
        settings.truncation_multiplier = self.confidence_radius
        settings.feasibility_cache_size = tfs_default_feasibility_cache_size
        settings.drg_epsilon = self.drg_epsilon
        settings.flux_epsilon = self.problem._epsilon_f[0][0]
        settings.min_rel_region_length = tfs_default_min_rel_region_length

        self.update_settings_function = lambda s, steps: tfs._fill_settings(
            s, self.num_samples, steps, points_array.shape[1], self.warmup, self.max_threads, self.num_dir_samples
        )
        logger.info("settings prepared: threads=%s, num_samples=%s", self.max_threads, self.num_samples)
        self.update_settings_function(settings, self.num_initial_steps)
        self.make_sampling_result_function = lambda r: tfs._make_sampling_result(self, r, self.num_samples)
        self.settings = settings

    def _conv_manager(self, directions_transform=None):
        # lazily set num_initial_steps default if not provided
        if self.num_initial_steps is None:
            self.num_initial_steps = (self.B.dimensionality ** 2 + self.num_samples)

        # Convergence manager depends on your ConvManager implementation
        self._convergence_manager = cm.ConvergenceManager(self, self.num_initial_steps, True, directions_transform)

    def sample(self, points_array):
        if not isinstance(points_array, np.ndarray):
            points_array = np.hstack(points_array)
        self.sampling_settings(points_array)
        self._conv_manager(directions_transform=self.directions_transform)

        result = self._convergence_manager.run(
            self.settings,
            self.update_settings_function,
            self.make_sampling_result_function,
            points_array,
            self.max_steps,
            default_max_psrf
        )
        maxpsrf = np.max(result.psrf)
        logger.info("Sampling converged: %s (max psrf=%.3f)", maxpsrf < 1.1, maxpsrf)
        return result
