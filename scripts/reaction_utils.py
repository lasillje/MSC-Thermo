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

def list_blocked_reactions(tmodel, condition: str, output_log: str, processes = 1):
    "Returns a list of blocked reactions. Does not remove the reactions from the model."

    blocked = find_blocked_reactions(tmodel, open_exchanges=False, processes = processes)

    write_to_log(output_log, f" - Found {len(blocked)} blocked reactions under {condition}")
    for rxn in blocked:
        write_to_log(output_log, f" --- Blocked reaction: {rxn}")

    print(blocked)
    return(blocked)

def list_and_remove_blocked_reactions(tmodel, condition: str, output_log: str, processes = 1):
    b = list_blocked_reactions(tmodel, condition, output_log, 1)
    tmodel.remove_reactions(b)
    for rxn in tmodel.reactions:
        thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)
    tmodel.update_thermo_info(fit_unknown_dfG0=True)
    return tmodel

def refine_subsystems(df):
    "Try to refine subsystems further based on BiGG database"

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