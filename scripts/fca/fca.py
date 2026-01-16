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

            model.remove_reactions([rxn], remove_orphans=True)

    for r in model.reactions:
        if len(r.metabolites) == 0:
            model.remove_reactions([r])

    print(f"Ended with {len(model.reactions)} reactions.")
    return model

def find_cofactor_loops(model):
    # Define groups of cofactors to "ignore" for the core stoichiometry
    cofactor_ids = {
        'nad_c', 'nadh_c', 'nadp_c', 'nadph_c', 
        'q8_c', 'q8h2_c', 'mqn8_c', 'mqn8h2_c', 
        'atp_c', 'adp_c', 'gtp_c', 'gdp_c', 'amp_c',
        'h_c', 'h_e', 'pi_c', 'h2o_c' # Also ignore protons/water
    }
    
    diamond_groups = {}

    for rxn in model.reactions:
        # Ignore boundaries/biomass/exchange/transport
        if rxn.boundary or "biomass" in rxn.id.lower() or len(rxn.metabolites) == 1 or len({m.compartment for m in rxn.metabolites}) > 1:
            continue
            
        core = {}
        for met, coeff in rxn.metabolites.items():
            if met.id not in cofactor_ids:
                core[met.id] = coeff
        
        if not core: continue # Skip reactions made ONLY of cofactors (like THD2)
        
        # Create a stable "signature"
        signature = tuple(sorted(core.items()))
        
        if signature not in diamond_groups:
            diamond_groups[signature] = []
        diamond_groups[signature].append(rxn.id)

    # Return only groups with multiple reactions (the diamonds)
    return {sig: ids for sig, ids in diamond_groups.items() if len(ids) > 1}

def collapse_all_loops(model, cofactor_loop_map):
    for signature, rxn_ids in cofactor_loop_map.items():
        # 1. Pick the first reaction as the 'Representative'
        rep_rxn = model.reactions.get_by_id(rxn_ids[0])
        others = rxn_ids[1:]
        
        print(f"Collapsing cofactor loop: {rxn_ids} -> Rep rxn: {rep_rxn.id}")
        
        for other_id in others:
            other_rxn = model.reactions.get_by_id(other_id)
            
            # 2. Transfer capacity (Bounds)
            # This ensures we don't reduce the maximum possible flux of the pathway
            rep_rxn.lower_bound += other_rxn.lower_bound
            rep_rxn.upper_bound += other_rxn.upper_bound


            if rep_rxn.upper_bound > 100:
                rep_rxn.upper_bound = 100
            if rep_rxn.lower_bound < -100:
                rep_rxn.lower_bound = -100
            
            # 3. Remove the redundant reaction
            model.remove_reactions([other_rxn])
            
        balance = rep_rxn.check_mass_balance()
        if balance: # Should be empty dict {}
            print(f"Warning: {rep_rxn.id} is unbalanced after merge: {balance}")
            
    return model    

def compress_model_full(model, remove_cofactor_loops = False):

    print("ONLY REACTIONS")

    if remove_cofactor_loops:
        cofactor_loops = find_cofactor_loops(model)
        for x in cofactor_loops:
            print(x)
        model = collapse_all_loops(model, cofactor_loops)

    rxn_map = find_coupled_reactions(model)
    model = collapse_coupled_reactions(model, rxn_map)

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