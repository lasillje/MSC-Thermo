"""Lightweight container for thermodynamic-space arrays.

ThermodynamicSpaceMod is a tiny dataclass-like container that mirrors the
attributes the rest of the package expects from the original `ThermodynamicSpace`
object. It intentionally does not perform heavy computations — that is handled
in the tbasis module.
"""

import numpy as np

class ThermodynamicSpaceMod:
    """Container for thermodynamic-space inputs.

    Parameters are plain numpy arrays (no pint units attached) to keep the
    module easy to test and fast to serialize.
    """

    def __init__(
        self,
        S_constraints: np.ndarray,
        reaction_ids,
        metabolites,
        drg0_prime_mean,
        drg0_prime_cov_sqrt,
        log_conc_cov,
        log_conc_mean,
    ):
        # store arrays and reshape vectors as column vectors where appropriate
        self.S_constraints = S_constraints
        self.reaction_ids = reaction_ids
        self.metabolites = metabolites
        self.drg0_prime_mean = np.asarray(drg0_prime_mean).reshape(-1, 1)
        self.drg0_prime_cov_sqrt = np.asarray(drg0_prime_cov_sqrt)
        self.log_conc_cov = np.asarray(log_conc_cov)
        self.log_conc_mean = np.asarray(log_conc_mean).reshape(-1, 1)
