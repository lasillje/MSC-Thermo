import os.path as path

from thermo_flux.core.model import ThermoModel
from thermo_flux.solver.gurobi import variability_analysis
from cobra.flux_analysis import find_blocked_reactions

from scripts.logger import write_to_log

from time import sleep

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import requests


from scripts.gen_model import constrain_bounds_fva

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
    tmodel.remove_reactions(b, remove_orphans=True)
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


def tfva_run_scenarios_one_model(tmodel, name, condition, output_folder, REMOVE_BLOCKED=False, APPLY_FVA=False, ONLY_WRITE=False, RUN_DIRECTLY=False, mipgap = 0.001, reaction_list=None):
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

    tmodel.m.Params.TimeLimit = 99999
    tmodel.m.Params.MIPGap = mipgap
    #tmodel.m.optimize()

    vars = []

    rxns = tmodel.reactions if reaction_list is None else reaction_list

    for rxn_id in rxns:

        rxn = tmodel.reactions.get_by_id(rxn_id.id)
        idx = tmodel.reactions.index(rxn)

        if RUN_DIRECTLY:

            tmodel.m = None  

            tmodel.objective = tmodel.reactions.biomass_EX  
            tmodel.add_TFBA_variables()

            tmodel.m.Params.TimeLimit = 99999
            tmodel.m.Params.MIPGap = mipgap

            v_var = tmodel.mvars["v"][0][idx]

            print(f"Running reaction: {rxn_id}, {idx}, {v_var}")
            gm = variability_analysis(tmodel, [v_var])

            if not ONLY_WRITE:
                gm.optimize()
                save_multiscenario_solutions(gm, output_folder, f"{name}_{rxn_id.id}_{idx}")
            else:
                gm.write(f"{output_folder}{path.sep}{idx}_{name}_{condition}_{rxn_id.id}.mps.gz")
        else:
           # v_var = tmodel.mvars["v"][0][idx]
            vars.append(tmodel.mvars["v"][0][idx])
            print(f"Added reaction: {rxn_id}, {idx}")
    print(f"Count: {len(vars)}")

    if not RUN_DIRECTLY:

        gm = variability_analysis(tmodel, vars)

        if not ONLY_WRITE:
            gm.optimize()
            save_multiscenario_solutions(gm, output_folder, name)
            gm.write(f"{output_folder}{path.sep}{name}_{condition}.sol")
        else:
            gm.write(f"{output_folder}{path.sep}{name}_{condition}.mps.gz")

def tfva_run_scenarios_one_model_mets(tmodel, name, condition, output_folder, OUTPUT_LOG, REMOVE_BLOCKED=True, APPLY_FVA=True, ONLY_WRITE=False):
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

    tmodel.m.Params.TimeLimit = 99999
    tmodel.m.Params.MIPGap = 0.001
    #tmodel.m.optimize()

    vars = []
    for i, met in enumerate(tmodel.metabolites):
        v = tmodel.mvars["ln_conc"][0][i]
        vars.append(v)
        print(f"Added {met.id} : { v }")

    gm = variability_analysis(tmodel, vars)
    if not ONLY_WRITE:
        gm.optimize()
    
        save_multiscenario_solutions(gm, output_folder, name)
        gm.write(f"{output_folder}{path.sep}{name}_{condition}_mets.sol")
    else:
        gm.write(f"{output_folder}{path.sep}{name}_{condition}_mets.mps.gz")

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
