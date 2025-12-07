"""Description of the space of thermodynamics-related quantities of a metabolic network.
"""

import logging
from typing import Dict, List, Optional, Set, Tuple, Union

import cobra
import numpy as np
import numpy.linalg as la
import scipy.linalg
from cobra.util.array import create_stoichiometric_matrix
from component_contribution.linalg import LINALG
from enkie import CompartmentParameters, Metabolite
from enkie.constants import R
from enkie.estimators import (EquilibratorGibbsEstimator,
                              GibbsEstimatorInterface)
from enkie.io.cobra import parse_metabolites
from enkie.utils import qrvector, qvector

from pta.commons import Q
from pta.concentrations_prior import ConcentrationsPrior
from pta.constants import default_min_eigenvalue_tds_basis
from pta.utils import (apply_transform, covariance_square_root,
                    get_candidate_thermodynamic_constraints, to_reactions_ids,
                    to_reactions_idxs)

from pta import ThermodynamicSpace

class ExtThermodynamicSpace(ThermodynamicSpace):
    """
    Construction, description and manipulation of the thermodynamic space of
    a metabolic network.

    Extended with the ability to provide custom drG0 data

    Parameters
    ----------
    S_constraints : np.ndarray
        Stoichiometric matrix of the reactions covered by thermodynamic constraints.
    reaction_ids : List[str]
        Identifiers of the reactions covered by thermodynamic constraints.
    metabolites : List[Metabolite]
        List describing the metabolites in the network.
    parameters : CompartmentParameters, optional
        The physiological parameters (pH, ionic strength, ...) of each compartment.
    concentrations : ConcentrationsPrior, optional
        Prior distributions for the metabolite concentrations.
    estimator : GibbsEstimatorInterface, optional
        Object used to estimate Gibbs free energies.
    dfg0_estimate : Optional[Tuple[Q, Q]], optional
        Estimate of formation energies (mean and a square root of the covariance matrix)
        in case the user wants to specify them manually. This is only used if
        :code:`estimator` is :code:`None`.
    """

    def __init__(
        self,
        S_constraints: np.ndarray,
        reaction_ids: List[str],
        metabolites: List[Metabolite],
        parameters: CompartmentParameters = None,
        concentrations: ConcentrationsPrior = None,
        estimator: GibbsEstimatorInterface = None,
        dfg0_estimate: Optional[Tuple[Q, Q]] = None,
        drg0_estimate: Optional[Tuple[Q, Q]] = None,
    ):
        self.ext_drg0_estimate = drg0_estimate
        if self.ext_drg0_estimate is not None:
            self.ext_drg0_mean = self.ext_drg0_estimate[0]
            self.ext_drg0_cov_sqrt = self.ext_drg0_estimate[1]

        super().__init__(S_constraints=S_constraints, reaction_ids=reaction_ids, metabolites=metabolites, parameters=parameters, concentrations=concentrations, estimator=estimator, dfg0_estimate=dfg0_estimate)

    @property
    def drg0_prime_mean(self) -> Q:
        """Gets the mean of the corrected standard reaction energies."""
        if self.ext_drg0_estimate is not None:
            return self.ext_drg0_mean 
        return self.S_constraints.T @ self.dfg0_prime_mean

    # @property
    # def drg0_prime_cov(self) -> Q:
    #     """Gets the covariance of the corrected standard reaction energies."""
    #     return self.drg0_prime_cov_sqrt @ self.drg0_prime_cov_sqrt.T

    @property
    def drg0_prime_cov_sqrt(self) -> Q:
        """Gets a square root of the covariance of the corrected standard
        reaction energies.
        """
        if self.ext_drg0_estimate is not None:
            return self.ext_drg0_cov_sqrt
        return self.S_constraints.T @ self.dfg0_prime_cov_sqrt

    @property
    def log_conc_mean(self):
        """Gets the mean of the log-concentrations."""
        return self._log_conc_mean

    @property
    def log_conc_cov(self):
        """Gets the covariance of the log-concentrations."""
        return self._log_conc_cov

    @property
    def dfg_prime_mean(self):
        """Gets the mean of the formation energies."""
        return self.dfg0_prime_mean + self.log_conc_mean * R * self.parameters.T()

    @property
    def dfg_prime_cov(self):
        """Gets the covariance of the formation energies."""
        return self.dfg0_prime_cov + self.log_conc_cov * (R * self.parameters.T()) ** 2

    @property
    def drg_prime_mean(self):
        """Gets the mean of the reaction energies."""
        return self.S_constraints.T @ self.dfg_prime_mean

    @property
    def drg_prime_cov(self):
        """Gets the covariance of the reaction energies."""
        return self.S_constraints.T @ self.dfg_prime_cov @ self.S_constraints
    
def needed_metabolites(model: cobra.Model, metabolites_namespace: str = None):
    """
    Returns a list of final metabolite IDs, use for aligning custom dfG0_* input data
    """
    constrained_rxns = get_candidate_thermodynamic_constraints(model, metabolites_namespace)

    constrained_rxn_ids = to_reactions_ids(constrained_rxns, model)
    #print(constrained_rxn_ids)
    constrained_rxn_idxs = to_reactions_idxs(constrained_rxns, model)

    S_constraints = create_stoichiometric_matrix(model)[:, constrained_rxn_idxs]
    is_metabolite_needed = np.sum(np.abs(S_constraints), axis=1) > 0
    needed_metabolites_ids = [
        model.metabolites[i].id for i in np.where(is_metabolite_needed)[0]
    ]

    return (constrained_rxn_ids, needed_metabolites_ids)