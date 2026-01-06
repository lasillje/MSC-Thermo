import cobra
import numpy as np
from scipy.linalg import null_space
import networkx as nx
import pandas as pd

def find_coupled_reactions(model):
    
    S = cobra.util.array.create_stoichiometric_matrix(model, array_type="DataFrame" )
    reaction_ids = S.columns.tolist()
    S_np = S.values
    N = null_space(S_np)

    pairs = {}

    tol = 1e-9

    for i in range(len(reaction_ids)):
        cur_reaction = reaction_ids[i]
        pairs.setdefault(cur_reaction, {})

        for j in range(i + 1, len(reaction_ids)):
            next_reaction = reaction_ids[j]
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
            
            priority_ids = ['biomass_EX']

            representative = sorted(list(rxn_set))[0]
            
            for p_id in priority_ids:
                if p_id in rxn_set:
                    representative = p_id
            
        
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