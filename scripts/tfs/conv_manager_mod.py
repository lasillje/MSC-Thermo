#-------------------------------------------------------------------------------
# Name:        conv_manager_mod.py

# Purpose:     modified convergence manager to set own direction transform for pta sampler
# 
# 
# Author:      Nathan Fargier
#
# Created:     07/11/2025

##version 1.0
#-------------------------------------------------------------------------------

"""Management of MCMC simulation to achieve desired convergence criteria."""
import logging
import time
from typing import Any, Callable

import _pta_python_binaries as pb
import numpy as np
from pta.sampling.commons import SamplerInterface, SamplingException, SamplingResult

logger = logging.getLogger(__name__)


class ConvergenceManager:
    """Class managing MCMC simulations to achieve desired convergence criteria.

    Parameters
    ----------
    sampler : Sampler
        Object implementing a specific MCMC sampler.
    num_initial_steps : int, optional
        Initial number of steps for simulations, by default -1 (automatic).
    samples_based_rounding : bool, optional
        If true, the direction sampling distribution will be adjusted after each
        iteration based on the distribution of samples so far, by default False
    """

    def __init__(
        self,
        sampler: SamplerInterface,
        num_initial_steps: int = -1,
        samples_based_rounding: bool = False,
        directions_transform=None ##transform for direction sampler
    ):
        self._sampler = sampler
        self._num_initial_steps = (
            num_initial_steps if num_initial_steps > 0 else sampler.dimensionality ** 2
        )
        self._round = samples_based_rounding
        self.directions_transform=np.identity(self._sampler.dimensionality) if directions_transform is None else directions_transform
    def run(
        self,
        settings: pb.SamplerSettings,
        update_settings_function: Callable[[pb.SamplerSettings, int], None],
        make_sampling_result_function: Callable[[Any], SamplingResult],
        initial_points: np.ndarray,
        max_steps: int,
        max_psrf: float,
    ) -> SamplingResult:
        """Run the sampler using this manager until the given PSRF or maximum number of
        steps is reached.

        Parameters
        ----------
        settings : pb.SamplerSettings
            The initial settings for the sampler.
        update_settings_function : Callable[[pb.SamplerSettings, int], None]
            Function for updating the settings with different numbers of steps.
        make_sampling_result_function: Callable[[Any], SamplingResult]
            Function for packing the result of the native sampler in a SamplingResult
            object.
        initial_points : np.ndarray
            The initial points for the simulation.
        max_steps : int
            Maximum number of steps to simulate.
        max_psrf : float
            Maximum PSRF to declare convergence.

        Returns
        -------
        SamplingResult
            The result of the sampler.

        Raises
        ------
        SamplingException
            If the native sampler fails.
        """
        assert (
            initial_points.shape[1] == settings.num_chains
        ), "Sampling requires the same number of initial points and chains."
        logger.info(
            f"Starting sampler with {settings.num_chains} chains and "
            f"{settings.max_threads} threads."
        )

        num_steps = self._num_initial_steps
        total_steps = 0
        iteration = 0
        converged = False
        directions_transform = self.directions_transform
        sampling_start_time = time.time()

        # Sample the space increasing the number of steps and optionally improving the
        # direction sampling distribution until one of the stopping conditions is
        # reached.
        while (not converged) and total_steps < max_steps:
            iteration += 1
            iteration_start_time = time.time()
            try:
                logger.info(
                    f"Running sampling iteration {iteration} with {num_steps} steps..."
                )
                update_settings_function(settings, num_steps)
                result = self._sampler.simulate(
                    settings, initial_points, directions_transform
                )
                total_steps += num_steps
                resulttmp = make_sampling_result_function(result)
                resulttmp.orthants.to_csv(f'tmp_save_orthants.csv')
            except Exception as e:
                logger.error(f"Sampling failed. {e}")
                raise SamplingException("Sampling failed.") from e
            else:
                iteration_time = time.time() - iteration_start_time
                step_rate = num_steps / iteration_time
                # Decide whether convergence was reached or we need to simulate the
                # chains for longer.
                psrf = self._sampler.compute_psrf(result)
                #
                if np.all(psrf <= max_psrf):
                    converged = True
                else:
                    if total_steps >= max_steps:
                        break ##to avoid svd if last iteration
                    # Update number of steps and initial points.
                    num_steps = min(2 * num_steps, max_steps - total_steps)
                    chains = self._sampler.get_chains(result)
                    initial_points = chains[-1, :, :]

                    # If needed update the directions transform as well.
                    if self._round:
                        _, S, Vh = np.linalg.svd(
                            np.hstack(chains).T, full_matrices=False
                        )
                        directions_transform = (
                            Vh.T
                            @ np.diag(S / np.min(S))
                            @ np.identity(self._sampler.dimensionality)
                        )

                logger.info(
                    f"Sampling iteration {iteration} completed in {iteration_time:.2f}"
                    f"s ({step_rate:.1f} steps/s). Max PSRF = {np.max(psrf):.2f}."
                )

        sampling_time = time.time() - sampling_start_time
        logger.info(
            f"Sampling completed in {sampling_time:.2f}s after a total of {total_steps}"
            f" steps. Convergence criteria were {'' if converged else 'not '}satisfied "
            f"(Max PSRF = {np.max(psrf):.3f}, PSRF threshold = {max_psrf:.3f})."
        )
        return make_sampling_result_function(result)
