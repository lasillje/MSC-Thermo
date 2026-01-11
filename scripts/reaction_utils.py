from pathlib import Path
import os.path as path
from cobra.flux_analysis import find_blocked_reactions
import cobra
import thermo_flux
import numpy as np
from thermo_flux.core.model import ThermoModel
from scripts.logger import write_to_log
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns  # optional but nicer
import scipy
import requests
import ast
from time import sleep
from thermo_flux.solver.gurobi import variability_analysis, variability_results, compute_IIS
import enkie
import pta
from typing import Any, Dict, Tuple, Union
from enkie.distributions import (
    LogNormalDistribution,
    LogUniformDistribution,
    distribution_from_string,
    distribution_to_string,
)
from enkie import CompartmentParameters
from pta import ConcentrationsPrior
import equilibrator_api
Q = equilibrator_api.Q_

from scripts.metabolite_utils import metabolite_to_bigg
from scipy.stats import norm
from scripts.gen_model import constrain_bounds_fva

from scripts.metabolite_utils import apply_met_tva

from scripts.fca.fca import find_coupled_reactions
from scripts.fca.fca import validate_coupling

def list_blocked_reactions(tmodel, condition: str, output_log: str, processes = 1, open_exch = False):
    "Returns a list of blocked reactions. Does not remove the reactions from the model."

    blocked = find_blocked_reactions(tmodel, open_exchanges = open_exch, processes = processes)

    write_to_log(output_log, f" - Found {len(blocked)} blocked reactions under {condition}")
    for rxn in blocked:
        write_to_log(output_log, f" --- Blocked reaction: {rxn}")

    print(blocked)
    print(len(blocked))
    return(blocked)

def list_and_remove_blocked_reactions(tmodel, condition: str, output_log: str, processes = 1):
    b = list_blocked_reactions(tmodel, condition, output_log, 1)
    tmodel.remove_reactions(b)
    for rxn in tmodel.reactions:
        thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
    tmodel.update_thermo_info(fit_unknown_dfG0=True)
    return tmodel

def refine_subsystems(df):
    "Try to refine subsystems further based on BiGG database "

    print("Refining subsystems...")

    API_MODEL = "iML1515"
    API_BASE = f"http://bigg.ucsd.edu/api/v2/models/{API_MODEL}/reactions/"

    df["RefinedSubsystem"] = None

    for i, rxn_id in enumerate(df["Abbrevation"]):
        rxn_id = str(rxn_id).strip()
        if not rxn_id:
            continue
        url = API_BASE + rxn_id
        try:
            print(f"Searching for refined subsystem for {rxn_id} at: {url}")
            r = requests.get(url, timeout=5)
            if r.status_code == 200:
                data = r.json()
                subsystem = data.get("subsystem")

                #if 'tRNA' in subsystem:
                #    print(f"TRNA: {i}, {rxn_id}")

                # Some outputs may have different formats apparently
                if not subsystem and "results" in data and len(data["results"]) > 0:
                    subsystem = data["results"][0].get("subsystem")

                if subsystem:
                    df.loc[i, "RefinedSubsystem"] = subsystem
                    print(f"[{i+1}] {rxn_id}: {subsystem}")
                else:
                    df.loc[i, "RefinedSubsystem"] = "Not found"
                    print(f"[{i+1}] {rxn_id}: subsystem missing")
            else:
                df.loc[i, "RefinedSubsystem"] = df.loc[i, "Subsystem"] #"Not in BiGG"
                print(f"[{i+1}] {rxn_id}: not found ({r.status_code})")

        except Exception as e:
            print(f"[{i+1}] {rxn_id}: error {e}")
            df.loc[i, "RefinedSubsystem"] = "Error"

        sleep(0.2)

def tfva_write_scenarios(tmodel, condition, output_folder, OUTPUT_LOG, lnc_unit="M", reactions_list = None, REMOVE_BLOCKED=True, APPLY_FVA=True):
    "Writes TFVA scenario files to the specified output folder, 1 file for each reaction"
    blocked_p = list_blocked_reactions(tmodel, condition, OUTPUT_LOG, 1, False)
    print(len(blocked_p))

    if REMOVE_BLOCKED:
        tmodel.remove_reactions(blocked_p, remove_orphans=True)

        if APPLY_FVA:

            print("Bounds before FVA: ")
            for x in tmodel.reactions:
                print(f"{x.lower_bound}, {x.upper_bound}")

            constrain_bounds_fva(tmodel, OUTPUT_LOG)

            print("Bounds after FVA: ")
            for x in tmodel.reactions:
                print(f"{x.lower_bound}, {x.upper_bound}")

        for rxn in tmodel.reactions:
            thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
        tmodel.update_thermo_info(fit_unknown_dfG0=True)
    else:
        for rxn in blocked_p:
            tmodel.reactions.get_by_id(rxn).lower_bound = 0
            tmodel.reactions.get_by_id(rxn).upper_bound = 0

    #apply_met_tva(tmodel, "hpc/WT-Glc_I_TFVA_Conc.mps.gz_objval.txt")

    tmodel.m = None  
    tmodel.objective = tmodel.reactions.biomass_EX  
    tmodel.add_TFBA_variables(lnc_unit=lnc_unit)

    tmodel.m.Params.TimeLimit = 20
    tmodel.m.optimize()

    folder = "Blocked_Removed" if REMOVE_BLOCKED else "Blocked_Restricted"
    suffix = "REMOVED" if REMOVE_BLOCKED else "RESTRICTED"
    suffix_fva = "FVA" if APPLY_FVA else ""

    rxn_ids = [r.id for r in tmodel.reactions]
    if reactions_list is not None:
        rxn_ids = reactions_list

    for rxn_id in rxn_ids:
        rxn = tmodel.reactions.get_by_id(rxn_id)
        idx = tmodel.reactions.index(rxn)

        v_var = [tmodel.mvars["v"][0][idx]]
        gm = variability_analysis(tmodel, v_var)
        gm.write(f"{output_folder}{path.sep}{condition}_{lnc_unit}_{rxn_id}_{idx}_{suffix}_{suffix_fva}_tfva.mps.gz")


def tfva_write_scenarios_one_model(tmodel, condition, output_folder, OUTPUT_LOG, lnc_unit="M", reactions_list = None, REMOVE_BLOCKED=True, APPLY_FVA=True):
    "Writes TFVA scenario files to the specified output folder, 1 file for each reaction"
    blocked_p = list_blocked_reactions(tmodel, condition, OUTPUT_LOG, 1, False)
    print(len(blocked_p))

    if REMOVE_BLOCKED:
        tmodel.remove_reactions(blocked_p, remove_orphans=True)

        if APPLY_FVA:

            print("Bounds before FVA: ")
            for x in tmodel.reactions:
                print(f"{x.lower_bound}, {x.upper_bound}")

            constrain_bounds_fva(tmodel, OUTPUT_LOG)

            print("Bounds after FVA: ")
            for x in tmodel.reactions:
                print(f"{x.lower_bound}, {x.upper_bound}")

        for rxn in tmodel.reactions:
            thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
        tmodel.update_thermo_info(fit_unknown_dfG0=True)
    else:
        for rxn in blocked_p:
            tmodel.reactions.get_by_id(rxn).lower_bound = 0
            tmodel.reactions.get_by_id(rxn).upper_bound = 0

    #apply_met_tva(tmodel, "hpc/WT-Glc_I_TFVA_Conc.mps.gz_objval.txt")

    tmodel.m = None  
    tmodel.objective = tmodel.reactions.biomass_EX  
    tmodel.add_TFBA_variables(lnc_unit=lnc_unit)

    tmodel.m.Params.TimeLimit = 20
    tmodel.m.optimize()

    folder = "Blocked_Removed" if REMOVE_BLOCKED else "Blocked_Restricted"
    suffix = "REMOVED" if REMOVE_BLOCKED else "RESTRICTED"
    suffix_fva = "FVA" if APPLY_FVA else ""

    rxn_ids = [r.id for r in tmodel.reactions]
    if reactions_list is not None:
        rxn_ids = reactions_list
    
    links = find_coupled_reactions(tmodel)

    dont_count = []
    for key in links:
        for x in links[key]:
            dont_count.append(x)
    x = set(dont_count)
    print(len(x))

    reactions = [r.id for r in tmodel.reactions if r.id not in dont_count]

    vars = []
    for rxn_id in reactions:
        rxn = tmodel.reactions.get_by_id(rxn_id)
        idx = tmodel.reactions.index(rxn)

        v_var = tmodel.mvars["v"][0][idx]
        vars.append(v_var)
        print(f"Added reaction: {rxn_id}, {idx}, {v_var}")
    print(f"Count: {len(vars)}")

    gm = variability_analysis(tmodel, vars)
    gm.write(f"{output_folder}{path.sep}{condition}_{lnc_unit}_{rxn_id}_{idx}_{suffix}_{suffix_fva}_tfva.mps.gz")


def save_multiscenario_solutions(m, output_folder, name):
    #save multiscenario solutions
    no_scenarios = m.NumScenarios
    if no_scenarios > 0:
        obj_val = {}
        obj_bound = {}
        optimal_bounds = {}
        MIPGaps = {}
        for i in range(0, no_scenarios, 2):
            rxn_idx = int(i / 2)
            # Minimization:
            m.params.ScenarioNumber = i
            m.update()
            ObjBound = m.ScenNObjBound
            ObjVal = m.ScenNObjVal
        #  print(rxn_idx, ObjBound, ObjVal)
            if ObjVal != 0:
                MIPGap = abs((ObjBound-ObjVal)/ObjVal)
            else:
                MIPGap = 0

            obj_val[rxn_idx] = [(-1) * ObjVal]
            obj_bound[rxn_idx] = [(-1) * ObjBound]
            MIPGaps[rxn_idx] = [MIPGap]

            if MIPGap <= m.params.MIPGap:
                optimal_bounds[rxn_idx] = [(-1) * ObjBound]
            else:
                optimal_bounds[rxn_idx] = [float('nan')]

            # Maximization:
            m.params.ScenarioNumber = i + 1
            m.update()
            ObjBound = m.ScenNObjBound
            ObjVal = m.ScenNObjVal
            if ObjVal != 0:
                MIPGap = abs((ObjBound-ObjVal)/ObjVal)
            else:
                MIPGap = 0
            obj_val[rxn_idx].append((+1) * m.ScenNObjVal)
            obj_bound[rxn_idx].append((1) * m.ScenNObjBound)
            MIPGaps[rxn_idx].append(MIPGap)
            if MIPGap <= 0.0001:
                optimal_bounds[rxn_idx].append((1) * ObjBound)
            else:
                optimal_bounds[rxn_idx].append(float('nan'))

        with open(f"{output_folder}/{name}_objval.txt", "w") as f:
            for k, val in obj_val.items():
                f.writelines(f"{k}: {val}\n")


def tfva_run_scenarios_one_model(tmodel, name, condition, output_folder, REMOVE_BLOCKED=False, APPLY_FVA=False):
    "Writes TFVA scenario files to the specified output folder, 1 file for each reaction"
    blocked_p = list_blocked_reactions(tmodel, condition, None, 1, False)
    print(len(blocked_p))

    if REMOVE_BLOCKED:
        tmodel.remove_reactions(blocked_p, remove_orphans=True)

        if APPLY_FVA:

            print("Bounds before FVA: ")
            for x in tmodel.reactions:
                print(f"{x.lower_bound}, {x.upper_bound}")

            constrain_bounds_fva(tmodel, None)

            print("Bounds after FVA: ")
            for x in tmodel.reactions:
                print(f"{x.lower_bound}, {x.upper_bound}")

        for rxn in tmodel.reactions:
            thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
        tmodel.update_thermo_info(fit_unknown_dfG0=True)
    else:
        for rxn in blocked_p:
            tmodel.reactions.get_by_id(rxn).lower_bound = 0
            tmodel.reactions.get_by_id(rxn).upper_bound = 0

    #apply_met_tva(tmodel, "hpc/WT-Glc_I_TFVA_Conc.mps.gz_objval.txt")

    tmodel.m = None  
    tmodel.objective = tmodel.reactions.biomass_EX  
    tmodel.add_TFBA_variables()

    tmodel.m.Params.TimeLimit = 500
    #tmodel.m.optimize()

    vars = []
    for rxn_id in tmodel.reactions:
        rxn = tmodel.reactions.get_by_id(rxn_id.id)
        idx = tmodel.reactions.index(rxn)

        v_var = tmodel.mvars["v"][0][idx]
        vars.append(v_var)
        print(f"Added reaction: {rxn_id}, {idx}, {v_var}")
    print(f"Count: {len(vars)}")

    gm = variability_analysis(tmodel, vars)

    gm.optimize()

    save_multiscenario_solutions(gm, output_folder, name)
    gm.write(f"{output_folder}{path.sep}{name}_{condition}.sol")

def tfva_run_scenarios_one_model_mets(tmodel, name, condition, output_folder, OUTPUT_LOG, REMOVE_BLOCKED=True, APPLY_FVA=True):
    "Writes TFVA scenario files to the specified output folder, 1 file for each reaction"
    blocked_p = list_blocked_reactions(tmodel, condition, OUTPUT_LOG, 1, False)
    print(len(blocked_p))

    if REMOVE_BLOCKED:
        tmodel.remove_reactions(blocked_p, remove_orphans=True)

        if APPLY_FVA:

            print("Bounds before FVA: ")
            for x in tmodel.reactions:
                print(f"{x.lower_bound}, {x.upper_bound}")

            constrain_bounds_fva(tmodel, OUTPUT_LOG)

            print("Bounds after FVA: ")
            for x in tmodel.reactions:
                print(f"{x.lower_bound}, {x.upper_bound}")

        for rxn in tmodel.reactions:
            thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
        tmodel.update_thermo_info(fit_unknown_dfG0=True)
    else:
        for rxn in blocked_p:
            tmodel.reactions.get_by_id(rxn).lower_bound = 0
            tmodel.reactions.get_by_id(rxn).upper_bound = 0

    tmodel.m = None  
    tmodel.objective = tmodel.reactions.biomass_EX  
    tmodel.add_TFBA_variables()

    tmodel.m.Params.TimeLimit = 500
    tmodel.m.optimize()

    vars = []
    for i, met in enumerate(tmodel.metabolites):
        v = tmodel.mvars["ln_conc"][0][i]
        vars.append(v)
        print(f"Added {met.id} : { v }")

    gm = variability_analysis(tmodel, vars)
    gm.optimize()
    
    save_multiscenario_solutions(gm, output_folder, name)
    gm.write(f"{output_folder}{path.sep}{name}_{condition}_mets.sol")

def tfva_update_bounds(tmodel, condition, tfva_results_dir):

    updated_rxns = set()

    directory = Path(tfva_results_dir)
    if not directory.exists():
        raise FileNotFoundError(f"Directory not found: {directory}")

    updated = 0
    skipped = 0
    not_found = 0

    for filepath in directory.glob("*.txt"):
        filename = filepath.name

        if condition and condition not in filename:
            continue

        parts = filename.split("_")
        if len(parts) < 5:
            print(f"Skipping incorrect filename: {filename}")
            skipped += 1
            continue

        rxn_id = parts[2]

        try:
            content = filepath.read_text().strip()
            # Remove unused prefix
            if ":" in content:
                bounds_str = content.split(":", 1)[1]
            else:
                bounds_str = content

            bounds = ast.literal_eval(bounds_str.strip())
            if not (isinstance(bounds, (list, tuple)) and len(bounds) == 2):
                raise ValueError("Not a [lb, ub] pair")
            lb, ub = float(bounds[0]), float(bounds[1])

        except Exception as e:
            print(f"Failed to parse bounds in {filename}: {e}")
            skipped += 1
            continue

        # Apply to model
        if rxn_id not in tmodel.reactions:
            print(f"Reaction not in model: {rxn_id} ({filename})")
            not_found += 1
            continue

        rxn = tmodel.reactions.get_by_id(rxn_id)
        old_lb, old_ub = rxn.lower_bound, rxn.upper_bound
        rxn.lower_bound = lb
        rxn.upper_bound = ub

        print(f"{rxn_id}: [{old_lb: .6f}, {old_ub: .6f}] -> [{lb: .6f}, {ub: .6f}]")
        updated += 1
        updated_rxns.add(rxn_id)

    print(f"\n=== TFVA Bounds Update Summary ===")
    print(f"Updated: {updated} reactions")
    print(f"Skipped/Malformed: {skipped}")
    print(f"Reaction not found in model: {not_found}")
    print(f"Total files processed: {updated + skipped + not_found}")
    print(f"Reactions in model which were not found in results: \n")
    missing_rxns = [x for x in tmodel.reactions if x not in updated_rxns]
    print(missing_rxns)

    return updated

def pta_get_parameters(tmodel, add_fake_p = False, temp_k = 310.15):
    
    # Our E. coli model has a temp of 310.15 (from model.xlsx)

    compartments = ["e", "c"]

    # Enkie does use a periplasm ('p') compartment which is missing in our model
    # Periplasm values will use the 'e' compartment's values with the same difference applied as the base enkie 'p' compartment
    # Default Enkie temp is 310.5 K
    # pH[p] = pH[e] + 0.01
    # pMg[p] = pMg[e]
    # I[p] = I[e]
    # phi[p] = phi[e] - 0.001
    # Also, we only have the difference of membrane potential (Phic - Phie), while Enkie requires both seperately
    # For now, phi[e] = 0, and phi[c] = tmodel.phi['ce'] * -1
    # To keep consistent with Enkie
    # https://gitlab.com/csb.ethz/enkie/-/blob/main/enkie/data/compartment_parameters/e_coli.csv?ref_type=heads
    
    phi_values = {"e": 0, "c": tmodel.phi['ce'].magnitude * -1}

    compartment_pH: Dict[str, Q] = {}
    compartment_pMg: Dict[str, Q] = {}
    compartment_I: Dict[str, Q] = {}
    compartment_phi: Dict[str, Q] = {}
    T = Q(temp_k, "K")

    for c in compartments:
        compartment_pH[c] = tmodel.pH[c]
        compartment_pMg[c] = tmodel.pMg[c]
        compartment_I[c] = tmodel.I[c]
        compartment_phi[c] = Q(phi_values[c], "V")
    
    if(add_fake_p):
        compartment_pH['p'] = compartment_pH['e'] + Q(0.01)
        compartment_pMg['p'] = compartment_pMg['e']
        compartment_I['p'] = compartment_I['e']
        compartment_phi['p'] = compartment_phi['e'] - Q(0.001, "V")

    return CompartmentParameters(
        compartment_pH, compartment_pMg, compartment_I, compartment_phi, T
    )

# def pta_get_concentrations(tmodel, metabolite_namespace = None,):
#     # ThermoModel metabolite concentrations are defined in linear space
#     # PTA requires them to be in ln space, so they will be converted to that
#     # Default data: https://gitlab.com/csb.ethz/pta/-/blob/main/pta/data/concentration_priors/M9_aerobic.csv?ref_type=heads

#     # Taken from pta concentrations_prior.py
#     default_log_conc = LogNormalDistribution(-8.3371, 1.9885)
#     #default_log_conc = LogNormalDistribution(-6.3371, 6.9885)

#     met_distributions = {}
#     default_distribution: Any = default_log_conc

#     # Create log normal distributions for every met
#     for met in tmodel.metabolites:
#         lb = met.lower_bound.magnitude
#         ub = met.upper_bound.magnitude

#         if(lb < ub or (lb == 0 and ub == 0)):
#             continue

#         met_id_split = met.id.split("_")
#         met_compartment = None

#         if(met_id_split):
#             for s in met_id_split:
#                 if s == "c" or s == "e":
#                     met_compartment = s

#         if met_compartment == None:
#             continue
        
#         log_dist = conc_to_logdist(lb, ub)
        
#         met_id_cleaned = metabolite_to_bigg(met.id)

#         met_name = met_id_cleaned if metabolite_namespace is None else f"{metabolite_namespace}:{met_id_cleaned}"
#         met_distributions[(met_name, met_compartment)] = log_dist

#     return ConcentrationsPrior(
#         metabolite_distributions = met_distributions,
#         default_distribution = default_distribution
#     )

def fix_negative_upper_bounds(model, tolerance: float = 1e-9):
    """
    Identifies boundary reactions that have a small, non-zero negative upper bound,
    which can cause a ValueError when internal PTA functions attempt to set 
    reaction.lower_bound = 0 before setting reaction.upper_bound = 0.
    
    If reaction.upper_bound is < 0, it means the reaction is forced to only
    operate in the reverse direction. To allow the boundary blocking step (v=0),
    we temporarily set the upper bound to 0

    An oversight in the PTA package causes it to error out when UB < 0, so this fix is necessary for now..
    """
    modified_count = 0
    print("\nStarting check for boundary reactions with negative upper bounds...")
    
    # A set of reactions where the original UB < 0
    reactions_to_fix = []
    
    for reaction in model.boundary:

        if reaction.upper_bound < -tolerance:
            reactions_to_fix.append(reaction)
            print(f"Found boundary reaction '{reaction.id}' with UB: {reaction.upper_bound:.4f}. Fixing...")
            
    for reaction in reactions_to_fix:
        # Reaction upper bounds are set to 0.0 to fix the mistake in PTA of setting the lower bound to 0 first, which causes a ValueError
        reaction.upper_bound = 0.0 
        modified_count += 1
    
    if modified_count == 0:
        print("No boundary reactions required fixing.")
    else:
        print(f"Completed boundary fix. Total reactions temporarily modified: {modified_count}")
        
    return reactions_to_fix

def count_blocked_pathways(reactions, name, condition, model_xlsx: str):
    "Counts and plots blocked pathways from a given list of blocked reactions and model xlsx file "
    df = pd.read_excel(model_xlsx, sheet_name="Reactions")

    refine_subsystems(df)

    df.columns = df.columns.str.strip()
    blocked_rxns = reactions

    df_blocked = df[df["Abbrevation"].isin(blocked_rxns)]
    subsystem_counts = df_blocked["RefinedSubsystem"].value_counts().reset_index()
    subsystem_counts.columns = ["RefinedSubsystem", "BlockedReactionCount"]

    print(subsystem_counts.head())

    #Plot absolute blocked
    plt.figure(figsize=(10, 6))
    sns.barplot(
        data=subsystem_counts,
        y="RefinedSubsystem",
        x="BlockedReactionCount",
        palette="viridis"
    )
    plt.title("Blocked Reactions per Subsystem")
    plt.xlabel("Number of Blocked Reactions")
    plt.ylabel("Subsystem")
    plt.tight_layout()
    plt.savefig(f"graphs{path.sep}blocked_abs_{name}_{condition}.png")
    plt.savefig(f"graphs{path.sep}blocked_abs_{name}_{condition}.svg")

    subsystem_total = df["RefinedSubsystem"].value_counts().reset_index()
    subsystem_total.columns = ["RefinedSubsystem", "TotalReactions"]

    merged = pd.merge(subsystem_total, subsystem_counts, on="RefinedSubsystem", how="left").fillna(0)
    merged["FractionBlocked"] = merged["BlockedReactionCount"] / merged["TotalReactions"]

    merged = merged.sort_values("FractionBlocked", ascending=False)

    #Plot fraction of blocked reactions
    plt.figure(figsize=(10, 6))
    sns.barplot(
        data=merged,
        y="RefinedSubsystem",
        x="FractionBlocked",
        palette="magma"
    )
    plt.xlabel("Fraction of Reactions Blocked")
    plt.ylabel("Subsystem")
    plt.title("Fraction of Reactions Blocked per Subsystem")
    plt.tight_layout()
    plt.savefig(f"graphs{path.sep}blocked_rel_{name}_{condition}.png")
    plt.savefig(f"graphs{path.sep}blocked_rel_{name}_{condition}.svg")

    return


def plot_calc_vs_exp(model_name: str, condition: str, sol_data: str, exp_data: str):

    gurobi_df = pd.read_csv(sol_data)
    exp_df = pd.read_csv(exp_data)

    gurobi_fluxes = gurobi_df[['reaction', 'v']].copy()

    exp_fluxes = exp_df[['cond', 'rxn', 'mean', 'sd']].copy()
    exp_fluxes = exp_fluxes[exp_fluxes['cond'] == condition]

    merged = pd.merge(exp_fluxes, gurobi_fluxes, left_on='rxn', right_on='reaction', how='inner')

    plt.figure(figsize=(7,6))
    ax = sns.scatterplot(
        data=merged,
        x='mean', 
        y='v', 
        hue='rxn',      # color by reaction
        s=50,          # dot size
        edgecolor='black'
    )

    plt.errorbar(
        merged['mean'], merged['v'], 
        xerr=merged['sd'], fmt='none', 
        ecolor='gray', alpha=0.6
    )

    ax.axhline(0, color='black', lw=1)
    ax.axvline(0, color='black', lw=1)
    plt.xlim(-20, 20)
    plt.ylim(-20, 20)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.xlabel("Measured flux (mean +/- sd)")
    plt.ylabel("Predicted flux")
    plt.title(f"Flux comparison for {condition}")
    plt.legend(title="Reaction", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

