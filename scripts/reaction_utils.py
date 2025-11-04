from cobra.flux_analysis import find_blocked_reactions
import cobra
import thermo_flux
from thermo_flux.core.model import ThermoModel
from scripts.logger import write_to_log

def list_blocked_reactions(tmodel, condition: str, output_log: str, processes = 1):
    "Returns a list of blocked reactions. Does not remove the reactions from the model."

    blocked = find_blocked_reactions(tmodel, processes = processes)
    write_to_log(output_log, f" - Found {len(blocked)} blocked reactions under {condition}")
    for rxn in blocked:
        write_to_log(output_log, f" --- Blocked reaction: {rxn}")
    return(blocked)