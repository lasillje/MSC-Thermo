from cobra import Model

def metabolite_to_bigg(met_id):
    """
    Takes a metabolite ID and changes it to a valid BiGG ID
    """
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