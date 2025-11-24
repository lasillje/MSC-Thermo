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
from time import sleep
from thermo_flux.solver.gurobi import variability_analysis, variability_results, compute_IIS

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

def tfva_write_scenarios(tmodel, condition, output_folder, OUTPUT_LOG, lnc_unit="M", num_reactions=-1, REMOVE_BLOCKED=True):
    "Writes TFVA scenario files to the specified output folder, 1 file for each reaction"
    blocked_p = list_blocked_reactions(tmodel, condition, OUTPUT_LOG, 1, False)
    print(len(blocked_p))

    if REMOVE_BLOCKED:
        tmodel.remove_reactions(blocked_p, remove_orphans=True)
        for rxn in tmodel.reactions:
            thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
        tmodel.update_thermo_info(fit_unknown_dfG0=True)
    else:
        for rxn in blocked_p:
            tmodel.reactions.get_by_id(rxn).lower_bound = 0
            tmodel.reactions.get_by_id(rxn).upper_bound = 0

    tmodel.m = None  
    tmodel.objective = tmodel.reactions.biomass_EX  
    tmodel.add_TFBA_variables(lnc_unit=lnc_unit)

    tmodel.m.Params.TimeLimit = 20
    tmodel.m.optimize()

    folder = "Blocked_Removed" if REMOVE_BLOCKED else "Blocked_Restricted"
    suffix = "B_REMOVED" if REMOVE_BLOCKED else "B_RESTRICTED"
    
    rxn_ids = [r.id for r in tmodel.reactions]
    max_rxns = len(rxn_ids) if num_reactions < 0 else num_reactions
    rxn_ids = rxn_ids[:max_rxns]

    for rxn_id in rxn_ids:
        rxn = tmodel.reactions.get_by_id(rxn_id)
        idx = tmodel.reactions.index(rxn)

        v_var = [tmodel.mvars["v"][0][idx]]
        gm = variability_analysis(tmodel, v_var)
        gm.write(f"{output_folder}{path.sep}{condition}_{lnc_unit}_{rxn_id}_{idx}_{suffix}_tfva.mps.gz")


def tfva_update_bounds(tmodel, condition, tfva_results_dir):
    print("wip")

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

