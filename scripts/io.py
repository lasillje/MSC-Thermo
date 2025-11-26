import cobra.io
import json
import numpy as np
from io import StringIO
from cobra import Model, Reaction, Metabolite
from scripts.metabolite_utils import metabolite_to_bigg

def _sanitize_value(value, replacement=0.0):
    if isinstance(value, float) and (np.isnan(value) or np.isinf(value)):
        return replacement
    return value

def export_tmodel_cobra(tmodel, output_path):
    """
    Exports a ThermoModel to a JSON file to load in later as a Cobra model
    """

    data = {
        "id": tmodel.id,
        "name": tmodel.name,
        "reactions": [],
        "metabolites": [],
        "genes": [g.id for g in tmodel.genes],
        "notes": tmodel.notes,
        "annotation": tmodel.annotation
    }

    for met in tmodel.metabolites:
        met_data = {
            "id": met.id,
            "name": met.name,
            "compartment": met.compartment,
            "charge": _sanitize_value(met.charge),
            "formula": met.formula,
        }
        data["metabolites"].append(met_data)


    #met_map = {met.id: met for met in tmodel.metabolites}

    print(data["metabolites"])

    for rxn in tmodel.reactions:
        lb_clean = _sanitize_value(rxn.lower_bound, replacement=-1000.0)
        ub_clean = _sanitize_value(rxn.upper_bound, replacement=1000.0)
        
        stoichiometry = {}
        for met_id, coeff in rxn.metabolites.items():
            stoichiometry[met_id.id] = _sanitize_value(coeff, replacement=0.0)

        rxn_data = {
            "id": rxn.id,
            "name": rxn.name,
            "lower_bound": lb_clean,
            "upper_bound": ub_clean,
            "gene_reaction_rule": rxn.gene_reaction_rule,
            "metabolites": stoichiometry,
            "subsystem": rxn.subsystem,
        }
        data["reactions"].append(rxn_data)

    with open(output_path, 'w') as f:
        json.dump(data, f, indent=4)
        
    print(f"Successfully exported clean COBRA data to {output_path}")
    return output_path


def load_tmodel_cobra(input_path, strip_compartment = True, metabolite_namespace = None):
    """
    Loads JSON file to create a cobra model
    """

    with open(input_path, 'r') as f:
        data = json.load(f)

    new_model = Model(data["id"])
    new_model.name = data["name"]
    new_model.notes = data.get("notes", {})
    new_model.annotation = data.get("annotation", {})

    metabolite_objects = {}
    for met_data in data["metabolites"]:
        met = Metabolite(**met_data)
        metabolite_objects[met.id] = met

    new_model.add_metabolites(list(metabolite_objects.values()))

    reaction_objects = []
    for rxn_data in data["reactions"]:
        rxn = Reaction(rxn_data["id"])
        rxn.name = rxn_data["name"]
        rxn.lower_bound = rxn_data["lower_bound"]
        rxn.upper_bound = rxn_data["upper_bound"]
        rxn.gene_reaction_rule = rxn_data["gene_reaction_rule"]
        rxn.subsystem = rxn_data["subsystem"]

        if rxn.lower_bound > rxn.upper_bound:
            print(f"Inconsistency fixed in {rxn.id}: LB ({rxn.lower_bound}) > UB ({rxn.upper_bound}). Clamping LB to UB.")
            rxn.lower_bound = rxn.upper_bound

        stoichiometry = {}
        for met_id, coeff in rxn_data["metabolites"].items():
            stoichiometry[metabolite_objects[met_id]] = coeff
        
        rxn.add_metabolites(stoichiometry)
        reaction_objects.append(rxn)

    new_model.add_reactions(reaction_objects)

    if metabolite_namespace is not None:
        for met in new_model.metabolites:
            met_name = met.id
            if strip_compartment:
                met_name = metabolite_to_bigg(met_name)
            met.annotation[metabolite_namespace] = met_name

    print(f"Successfully loaded clean COBRA model: {new_model.id}")
    return new_model