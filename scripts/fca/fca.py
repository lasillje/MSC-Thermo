import cobra
import numpy as np
from scipy.linalg import null_space
import networkx as nx
import pandas as pd
from equilibrator_api import Q_

def find_coupled_reactions(model):
    
    protected_ids = set()
    for rxn in model.reactions:
        if "biomass" in rxn.id.lower():
            protected_ids.add(rxn.id)
        if rxn.boundary or len(rxn.metabolites) == 1:
            protected_ids.add(rxn.id)
        if len({m.compartment for m in rxn.metabolites}) > 1:
            protected_ids.add(rxn.id)

    S = cobra.util.array.create_stoichiometric_matrix(model, array_type="DataFrame" )
    reaction_ids = S.columns.tolist()
    S_np = S.values
    N = null_space(S_np)

    pairs = {}

    tol = 1e-9

    for i in range(len(reaction_ids)):
        cur_reaction = reaction_ids[i]
        if cur_reaction in protected_ids: continue
        pairs.setdefault(cur_reaction, {})

        for j in range(i + 1, len(reaction_ids)):
            next_reaction = reaction_ids[j]
            if next_reaction in protected_ids: continue
            pairs.setdefault(next_reaction, {})

            n_i = N[i, :]
            n_j = N[j, :]

            active_columns = np.logical_and(np.abs(n_i) > tol, np.abs(n_j) > tol)

            if np.any(active_columns):
                ratios = n_i[active_columns] / n_j[active_columns]
                
                if np.max(ratios) - np.min(ratios) < tol:
                    ratio = ratios[0]

                    pairs[cur_reaction][next_reaction] = ratio
                    pairs[next_reaction][cur_reaction] = 1 / ratio


    G = nx.Graph()
    G.add_nodes_from(reaction_ids)

    for rxn_a, targets in pairs.items():
        for rxn_b in targets.keys():
            if rxn_a < rxn_b:
                G.add_edge(rxn_a, rxn_b)
    pos = nx.spring_layout(G, k=0.15, iterations=5000)
    nx.draw_networkx_edges(G, pos)
    coupled_reaction_sets = list(nx.connected_components(G))

    final_rxns = {}

    for rxn_set in coupled_reaction_sets:
        if len(rxn_set) > 1:
            
            representative = sorted(list(rxn_set))[0]
        
            linked_reactions = {}
            for rxn in rxn_set:
                if rxn == representative:
                    continue
                
                rep_index = reaction_ids.index(representative)
                rxn_index = reaction_ids.index(rxn)
                
                max_rep_index = np.argmax(np.abs(N[rep_index, :]))
                
                ratio_rxn_to_rep = N[rxn_index, max_rep_index] / N[rep_index, max_rep_index]

                linked_reactions[rxn] = ratio_rxn_to_rep
            
            final_rxns[representative] = linked_reactions
            
    return final_rxns

def find_dead_end_metabolites(model):
    """
    Identifies metabolites that are only produced or only consumed.
    """
    dead_ends = []
    for metabolite in model.metabolites:
        rxns = metabolite.reactions
        
        if len(rxns) == 0:
            dead_ends.append(metabolite)
            continue
            
        #coeffs = [r.get_coefficient(metabolite) for r in rxns]
        
        #if all(c > 0 for c in coeffs) or all(c < 0 for c in coeffs):
        #    dead_ends.append(metabolite)
    
    print(f"Found {len(dead_ends)} dead-end metabolites.")
    for x in dead_ends:
        print(x.id)
    return dead_ends

def find_coupled_metabolites(model):
    
    tol = 1e-9

    S = cobra.util.array.create_stoichiometric_matrix(model, array_type="DataFrame" )

    proportional_rows = {}

    for met_id in S.index:

        row = S.loc[met_id].values

        non_zeros = np.abs(row) > tol

        if not np.any(non_zeros):
            continue

        indices = np.where(non_zeros)[0]
        values = row[indices]

        factor = values[0]
        normalized_values = values / factor

        row_hash = (tuple(indices), tuple(np.round(normalized_values, 6)))

        if row_hash not in proportional_rows:
            proportional_rows[row_hash] = []

        proportional_rows[row_hash].append({'id': met_id, 'factor': factor})

    coupling_map = {}

    for hash, group in proportional_rows.items():
        group_sorted = sorted(group, key=lambda x: x['id'])

        representative = group_sorted[0]
        
        rep_id = representative['id']
        rep_factor = representative['factor']

        linked_mets = {}

        for item in group_sorted[1:]:
            ratio = item['factor'] / rep_factor
            linked_mets[item['id']] = ratio

        coupling_map[rep_id] = linked_mets

    print(coupling_map)
    return coupling_map

def collapse_coupled_reactions(model, coupling_map):
    tol = 1e-9
    print(f"Excluding biomass, boundary reactions.")
    print(f"Starting with {len(model.reactions)} reactions.")
    for rep, linked in coupling_map.items():
        rep_rxn = model.reactions.get_by_id(rep)

        for rxn_id, ratio in linked.items():
            rxn = model.reactions.get_by_id(rxn_id)

            target_metabolites = {met: coeff * ratio for met, coeff in rxn.metabolites.items()}
            rep_rxn.add_metabolites(target_metabolites)

            for met in list(rep_rxn.metabolites):
                if abs(rep_rxn.get_coefficient(met)) < tol:
                    rep_rxn.add_metabolites({met: -rep_rxn.get_coefficient(met)})

            if ratio > 0:
                rep_rxn.lower_bound = max(
                    rep_rxn.lower_bound, rxn.lower_bound / ratio
                )
                rep_rxn.upper_bound = min(
                    rep_rxn.upper_bound, rxn.upper_bound / ratio
                )
            else:
                rep_rxn.lower_bound = max(
                    rep_rxn.lower_bound, rxn.upper_bound / ratio
                )
                rep_rxn.upper_bound = min(
                    rep_rxn.upper_bound, rxn.lower_bound / ratio
                )

            model.remove_reactions([rxn])

    for r in model.reactions:
        if len(r.metabolites) == 0:
            model.remove_reactions([r])

    print(f"Ended with {len(model.reactions)} reactions.")
    return model

def collapse_coupled_metabolites(model, coupling_map):

    print(f"Starting with {len(model.metabolites)} metabolites.")

    dfg_dict = {met.id: getattr(met, 'dfG0prime') for met in model.metabolites}
    new_dfg_dict = dict()

    for rep_id, linked in coupling_map.items():

        base_energy = float(dfg_dict[rep_id].m)
        added_energy = 0.0

        for linked_id, ratio in linked.items():
            added_energy += ratio * float(dfg_dict[linked_id].m)
        
        new_dfg_dict[rep_id] = Q_(base_energy + added_energy, "kJ/M")
        model.metabolites.get_by_id(rep_id).dfG0prime = new_dfg_dict[rep_id]

    mets_to_remove = [model.metabolites.get_by_id(mid) for mid in linked.keys()]
    model.remove_metabolites(mets_to_remove)

    return model
    
def compress_model_full(model):

    print("ONLY REACTIONS")
    rxn_map = find_coupled_reactions(model)
    model = collapse_coupled_reactions(model, rxn_map)

    # If enabled dont forget to assign dfG0 before this point
    #met_map = find_coupled_metabolites(model)
    #model = collapse_coupled_metabolites(model, met_map)

    #dead_mets = find_dead_end_metabolites(model)
    #model.remove_metabolites(dead_mets)

    #print(f"Ended with {len(model.metabolites)} metabolites.")

    return model

def validate_coupling(cobra_model, coupling_results, samples=500, tolerance=1e-9):
    cobra_model.optimize()
    success = 0
    print(len(cobra_model.reactions))
    try:
        sampling_solution = cobra.sampling.sample(cobra_model, n=samples, processes=1)
    except Exception as e:
        print("Model is infeasible")
        return pd.DataFrame()
    
    validation_data = []

    # Iterate through all coupled pairs
    for rep_id in coupling_results:

        links = coupling_results[rep_id]
        print(f"Current: {rep_id}, {links}")
        if rep_id not in sampling_solution.columns:
            print(f"Reporting reaction {rep_id} not found in model flux space. Skipping.")
            continue

        rep_fluxes = sampling_solution[rep_id]

        for linked_id in links:
            expected_ratio = coupling_results[rep_id][linked_id]
            print(f"{linked_id} - {expected_ratio}")
            if linked_id not in sampling_solution.columns:
                print(f"Linked reaction {linked_id} not found. Skipping.")
                continue

            linked_fluxes = sampling_solution[linked_id]


            valid_indices = np.abs(rep_fluxes) > tolerance
            
            if not np.any(valid_indices):
                validation_data.append({
                    'Reporting_Rxn': rep_id,
                    'Linked_Rxn': linked_id,
                    'Expected_Ratio': expected_ratio,
                    'Validity': 'Skipped (flux always zero)',
                    'Pass_Rate': 0.0,
                    'Max_Error': np.nan
                })
                continue
            
            valid_rep_fluxes = rep_fluxes[valid_indices]
            valid_linked_fluxes = linked_fluxes[valid_indices]

            observed_ratios = valid_linked_fluxes / valid_rep_fluxes
            
            absolute_errors = np.abs(observed_ratios - expected_ratio)
            
            pass_count = np.sum(absolute_errors < tolerance)
            total_valid_samples = len(valid_rep_fluxes)
            
            pass_rate = (pass_count / total_valid_samples) if total_valid_samples > 0 else 0.0
            
            # Determine overall validity
            if pass_rate > 0.99:
                validity = 'Validated'
                success += 1
            elif pass_rate > 0.9:
                validity = 'Weak (Check Tolerance)'
            else:
                validity = 'Failed'

            validation_data.append({
                'Reporting_Rxn': rep_id,
                'Linked_Rxn': linked_id,
                'Expected_Ratio': expected_ratio,
                'Validity': validity,
                'Pass_Rate': pass_rate,
                'Max_Error': absolute_errors.max()
            })
    print(success)
    print(f"Passed: {success} / {len(cobra_model.reactions)}")
    return pd.DataFrame(validation_data)