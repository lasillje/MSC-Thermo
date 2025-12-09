from cobra import Model
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt
import pta
from pta import ConcentrationsPrior
import ast
from equilibrator_api import Q_

def metabolite_to_bigg(met_id):
    """
    Takes a metabolite ID and changes it to a valid BiGG ID
    """

    if "co2tot" in met_id:
        return "hco3"

    met_id_cleaned = met_id
    if "_c" in met_id or "_e" in met_id:
         met_id_cleaned = met_id_cleaned[:-2]
    
    # TODO: Make this better instead of a wall of if statements
    # Fix for amino acids
    if "-" in met_id_cleaned and met_id_cleaned[-2] == "-":
        if met_id_cleaned[-1] == "L" or met_id_cleaned[-1] == "D" or met_id_cleaned[-1] == "S" or met_id_cleaned[-1] == "M" or met_id_cleaned[-1] == "R":
            met_id_cleaned = met_id_cleaned.replace("-", "__")
        elif met_id_cleaned[-1] == "B" or met_id_cleaned[-1] == "T" or met_id_cleaned[-1].isnumeric():
            # Special case for trans/beta which have prefix _B and _T such as acon_T, dad_2
            met_id_cleaned = met_id_cleaned.replace("-", "_")

    # Special case for adphep_DD, 26dap_LL, etc.
    if "-" in met_id_cleaned and met_id_cleaned[-3] == "-":
        if met_id_cleaned[-2:] == "LD" or met_id_cleaned[-2:] == "DD" or met_id_cleaned[-2:] == "LL" or met_id_cleaned[-2:] == "DL":
            met_id_cleaned = met_id_cleaned.replace("-", "_")

    # Special case for peptidoEC, pcEC, lpsEC which are listed in BiGG as pc_EC, etc.
    if "EC" in met_id_cleaned and met_id_cleaned.endswith("EC"):
        met_id_cleaned = met_id_cleaned[:-2] + '_' + met_id_cleaned[-2:]

    return met_id_cleaned

def conc_to_logdist(lower_bound, upper_bound, frac_of_area = 0.95):
    
    if lower_bound <= 0 or upper_bound <= 0:
        raise ValueError("Concentration bounds must be positive")
    if  lower_bound > upper_bound:
        raise ValueError("Lower bound cannot be higher than upper bound")
    if frac_of_area <= 0.0 or frac_of_area > 1.0:
        raise ValueError("Fraction of area must be between 0.0 and 1.0")
    
    # Linear to ln
    ln_lower = np.log(lower_bound)
    ln_upper = np.log(upper_bound)

    # Standard deviations based on frac_of_area
    beta = (1.0 - frac_of_area) / 2.0 # Area per tail
    Z = norm.ppf(1.0 - beta) # About 1.96 for frac = 0.95

    log_mean = (ln_lower + ln_upper) /  2.0
    log_std = (ln_upper - ln_lower) / (2.0 * Z) # Distance of bounds divided by total number of standard deviations to get a std value

    log_dist = pta.LogNormalDistribution(log_mean=log_mean, log_std=log_std)

    return log_dist

def remove_orphan_metabolites(model: Model):
    """
    Removes metabolites from the model that do not participate in any reaction
    """

    linked_metabolites = set()
    for rxn in model.reactions:
        for met in rxn.metabolites:
            linked_metabolites.add(met)
            
    all_metabolites = set(model.metabolites)
    orphan_metabolites = list(all_metabolites - linked_metabolites)
    
    orphan_ids = [m.id for m in orphan_metabolites]
    
    if orphan_metabolites:
        print(f"\nRemoving {len(orphan_metabolites)} orphan metabolites:")
        print(f"Orphaned IDs: {', '.join(orphan_ids)}")
        
        model.remove_metabolites(orphan_metabolites)
        
    else:
        print("\nNo orphan metabolites found. ")
        
    return len(orphan_metabolites)

def apply_met_tva(tmodel, met_tva_file):
    bounds_dict = dict()
    with open(met_tva_file, "r") as f:
        for line in f:
            clean_line = line.strip()
            if not clean_line:
                print(f"Skipping line {clean_line}")
                continue
            try:
                index_str, bounds_str = clean_line.split(':', 1)
                index = int(index_str.strip())

                cleaned_bounds_str = bounds_str.strip().strip('[] ')
                lower_str, upper_str = cleaned_bounds_str.split(',')

                lower = float(lower_str.strip())
                upper = float(upper_str.strip())

                linear_lower = np.exp(lower) * 1e3 # From molar conc back to millimolar
                linear_upper = np.exp(upper) * 1e3 # Same
                
                bounds_dict[index] = [linear_lower, linear_upper]
            except ValueError as e:
                print(f"Skipping line due to parsing error: '{line.strip()}' - Error: {e}")
            except Exception as e:
                print(f"An unexpected error occurred while processing line: '{line.strip()}' - Error: {e}") 
    
    for met in tmodel.metabolites:
        met_index = tmodel.metabolites.index(met)

        if met_index not in bounds_dict:
            print(f"Skipped metabolite {met.id} as it was not found in TVA data.")
            continue

        cur_lower, cur_upper = met.lower_bound, met.upper_bound
        new_lower, new_upper = bounds_dict[met_index][0], bounds_dict[met_index][1]

        print(f"Metabolite {met.id} - Old: {cur_lower}, {cur_upper} | New: {new_lower, new_upper}")
    
        met.upper_bound = Q_(new_upper, "millimolar")
        met.lower_bound = Q_(new_lower, "millimolar")


def graph_ln_dist(mu, sigma):

    C_min_target = 1e-9
    C_max_target = 1e3
    ln_C_min_target = np.log(C_min_target)
    ln_C_max_target = np.log(C_max_target)

    ln_samples = norm.rvs(loc=mu, scale=sigma, size=1000)
    concentration_samples = np.exp(ln_samples)

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f'LogNormal Concentration Prior (mean={mu}, std={sigma})', fontsize=16)

    axes[0].hist(ln_samples, bins=50, density=True, color='skyblue', edgecolor='black', alpha=0.7)
    axes[0].axvline(mu, color='red', linestyle='--', label=f'Mean (μ) = {mu:.3f}')
    axes[0].axvline(ln_C_min_target, color='green', linestyle=':', label=f'ln(C_min) = {ln_C_min_target:.2f}')
    axes[0].axvline(ln_C_max_target, color='green', linestyle=':', label=f'ln(C_max) = {ln_C_max_target:.2f}')
    axes[0].set_title('Distribution of Log Concentration (ln(C))')
    axes[0].set_xlabel('ln(Concentration) [ln(M)]')
    axes[0].set_ylabel('Probability Density')
    axes[0].legend()
    axes[0].grid(True, linestyle='--')

    axes[1].hist(concentration_samples, bins=np.logspace(np.log10(C_min_target / 10), np.log10(C_max_target * 10), 50), 
             density=True, color='lightcoral', edgecolor='black', alpha=0.7)
    axes[1].axvline(C_min_target, color='green', linestyle=':', label=f'C_min = {C_min_target:.0e}')
    axes[1].axvline(C_max_target, color='green', linestyle=':', label=f'C_max = {C_max_target:.0e}')
    axes[1].set_xscale('log')
    axes[1].set_title('Distribution of Concentration (C)')
    axes[1].set_xlabel('Concentration [M] (Log Scale)')
    axes[1].set_ylabel('Probability Density')
    axes[1].legend()
    axes[1].grid(True, linestyle='--')

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()