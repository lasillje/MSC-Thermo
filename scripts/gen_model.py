
import os.path as path
import glob, os
from datetime import datetime
from importlib.metadata import version
import cobra
import thermo_flux
from thermo_flux.io import load_excel as ex
from thermo_flux.core.model import ThermoModel
from equilibrator_api import  Q_
import pandas as pd
from thermo_flux.io import helper_load as hl
import numpy as np
from thermo_flux.io import load_excel as ex
import gurobipy as gp
from scripts.logger import write_to_log
from gurobipy import GRB

def gen_model(name: str, model_xlsx: str, kegg: str, reed: str, inchi:str, gams: str, output_log: str, add_o2: bool, add_co2: bool):
    "Generates a base Thermo-Flux model given the input files"
    tmodel = ex.create_thermo_model(name, model_excel=model_xlsx, keggids_csv=kegg, edit_mets={})
    write_to_log(output_log, f"Loaded model {name} from {model_xlsx}")

    # ADD HYDROXYBENZOATE TRANSPORT & EXCHANGE:

    # Define extracellular orotate:
    hbz_e = cobra.Metabolite(id="4hbz_e", compartment="e")
    hbz_e = thermo_flux.core.metabolite.ThermoMetabolite(hbz_e, model=tmodel)
    hbz_e.annotation = tmodel.metabolites.get_by_id("4hbz_c").annotation  
    tmodel.metabolites.append(hbz_e)

    # Hydroxybenzoate H+ antiporter (aaeB):
    HBZt3 = cobra.Reaction("HBZt3")
    HBZt3.lower_bound = -1000
    HBZt3.upper_bound = 1000
    HBZt3.add_metabolites({tmodel.metabolites.get_by_id("4hbz_e"): +1,
                      tmodel.metabolites.h_e: -1,
                      tmodel.metabolites.get_by_id("4hbz_c"): -1,
                      tmodel.metabolites.h_c: +1})
    HBZt3 = thermo_flux.core.reaction.ThermoReaction(HBZt3, model=tmodel)

    # Exchange:
    EX_4hbz = cobra.Reaction("EX_4hbz")
    EX_4hbz.add_metabolites({tmodel.metabolites.get_by_id("4hbz_e"): -1})
    EX_4hbz = thermo_flux.core.reaction.ThermoReaction(EX_4hbz, model=tmodel)

    # Add reactions to the model
    tmodel.add_reactions([HBZt3, EX_4hbz])

    write_to_log(output_log, f"Added reactions for hydroxybenzoate transport and exchange")
    write_to_log(output_log, f"HBZT3    reaction: {tmodel.reactions.HBZt3}")
    write_to_log(output_log, f"EX_4hbz  reaction: {tmodel.reactions.HBZt3}")

    # ADD ISOPROPYLMALATE TRANSPORT & EXCHANGE:

    # Define extracellular orotate:
    ipm_e = cobra.Metabolite(id="3c3hmp_e", compartment="e")
    ipm_e = thermo_flux.core.metabolite.ThermoMetabolite(ipm_e, model=tmodel)
    ipm_e.annotation = tmodel.metabolites.get_by_id("3c3hmp_c").annotation  
    tmodel.metabolites.append(ipm_e)

    # Diffusion across the membrane:
    ipm_diff = cobra.Reaction("IPMex")
    ipm_diff.lower_bound = -1000
    ipm_diff.upper_bound = 1000
    ipm_diff.add_metabolites({tmodel.metabolites.get_by_id("3c3hmp_c"): -1,
                            tmodel.metabolites.get_by_id("3c3hmp_e"): +1})
    ipm_diff = thermo_flux.core.reaction.ThermoReaction(ipm_diff, model=tmodel)


    # Exchange:
    EX_3c3hmp = cobra.Reaction("EX_3c3hmp")
    EX_3c3hmp.add_metabolites({tmodel.metabolites.get_by_id("3c3hmp_e"): -1})
    EX_3c3hmp = thermo_flux.core.reaction.ThermoReaction(EX_3c3hmp, model=tmodel)

    # Add reactions to the model
    tmodel.add_reactions([ipm_diff, EX_3c3hmp])

    write_to_log(output_log, f"Added reactions for isopropyl-malate transport and exchange")
    write_to_log(output_log, f"IPMex     reaction: {tmodel.reactions.IPMex}")
    write_to_log(output_log, f"EX_3c3hmp reaction: {tmodel.reactions.EX_3c3hmp}")

    # ADD OROTATE TRANSPORT AND EXCHANGE:

    # Define extracellular orotate:
    oro_e = cobra.Metabolite(id="orot_e", compartment="e")
    oro_e = thermo_flux.core.metabolite.ThermoMetabolite(oro_e, model=tmodel)
    oro_e.annotation = tmodel.metabolites.orot_c.annotation  
    tmodel.metabolites.append(oro_e)

    # Diffusion across the membrane:
    oro_diff = cobra.Reaction("OROTex")
    oro_diff.lower_bound = -1000
    oro_diff.upper_bound = 1000
    oro_diff.add_metabolites({tmodel.metabolites.orot_c: -1,
                            tmodel.metabolites.orot_e: +1})
    oro_diff = thermo_flux.core.reaction.ThermoReaction(oro_diff, model=tmodel)

    # Dicarboxylate/H+ symporter (dctA):
    dcta = cobra.Reaction("DCTA")
    dcta.lower_bound = -1000
    dcta.upper_bound = 0
    dcta.add_metabolites({tmodel.metabolites.orot_e: +1,
                        tmodel.metabolites.h_e: +1,
                        tmodel.metabolites.orot_c: -1,
                        tmodel.metabolites.h_c: -1})
    dcta = thermo_flux.core.reaction.ThermoReaction(dcta, model=tmodel)

    # Exchange:
    EX_oro = cobra.Reaction("EX_oro")
    EX_oro.add_metabolites({tmodel.metabolites.orot_e: -1})
    EX_oro = thermo_flux.core.reaction.ThermoReaction(EX_oro, model=tmodel)

    # Add reactions to the model
    tmodel.add_reactions([oro_diff, dcta, EX_oro])

    write_to_log(output_log, f"Added reactions for orotate transport and exchange")
    write_to_log(output_log, f"OROTex     reaction: {tmodel.reactions.OROTex}")
    write_to_log(output_log, f"DCTA       reaction: {tmodel.reactions.DCTA}")
    write_to_log(output_log, f"EX_oro     reaction: {tmodel.reactions.EX_oro}")

    tmodel.pH = {"c": Q_(7.6), "e": Q_(7)} #pH
    tmodel.I = {"c": Q_(0.25,'M'), "e": Q_(0.25,'M')} #ionic stength
    tmodel.phi = {'ce':Q_(0.15,'V')} #membrane potential ‘ce’ represents the voltage between compartment ‘c’ and compartment 'e’ defined as Phic - Phie
    tmodel.pMg = {'e': Q_(3), 'c': Q_(3)}

    # Update metabolite annotations with the IDs from KEGG
    for met in tmodel.metabolites:
        met.annotation["bigg.metabolite"] = met.id[:-2]
    write_to_log(output_log, f"Updated metabolite annotations with KEGG IDs")

    # Update the inchi strings of some unknown metabolites
    df = pd.read_csv(reed, header=None).set_index(0)
    write_to_log(output_log, f"Reading InChi strings from: {reed}")

    unknown_mets = []
    for met in tmodel.metabolites:
        if met.id[:-2] in df.index:
            met.annotation['InChI'] = df.loc[met.id[:-2]].iloc[1]
            unknown_mets.append(met)
    write_to_log(output_log, f" - updated InChi strings (1/2)")

    for met in tmodel.metabolites:
        if ('kegg' in met.annotation):
            if (met.annotation['kegg'] in df.index):
            
                inchi = df.loc[met.annotation['kegg']].iloc[0]

                if type(inchi) is str:
                    met.annotation['InChI'] = df.loc[met.annotation['kegg']].iloc[0]
                    unknown_mets.append(met)
    write_to_log(output_log, f" - updated InChi strings (2/2)")

    # Additional data from excel spredsheet needs to be imported for some compounds with unknown structure
    sheets = ex.read_excelfile(model_xlsx)
    write_to_log(output_log, f"Reading model sheet from: {model_xlsx}")

    for met in tmodel.metabolites:
        if met.id[:-2] in sheets['Metabolites']['Unnamed: 0'].values:

            common_name = sheets['Metabolites'].loc[sheets['Metabolites']['Unnamed: 0'] == met.id[:-2]]['Unnamed: 11'].values

            protons = sheets['Metabolites'].loc[sheets['Metabolites']['Unnamed: 0'] == met.id[:-2]]['Unnamed: 9'].values

            charge =  sheets['Metabolites'].loc[sheets['Metabolites']['Unnamed: 0'] == met.id[:-2]]['Unnamed: 8'].values

            formula =  sheets['Metabolites'].loc[sheets['Metabolites']['Unnamed: 0'] == met.id[:-2]]['Unnamed: 12'].values

    #         print(met.id, protons[0])
            met.formula = 'H'+str(protons[0])

            if len(common_name) > 0:
                met.notes['common name'] = common_name[0]
                met.charge = charge[0]
            else:
                met.notes['common name'] = ''
    write_to_log(output_log, f"Updated metabolite annotations with charge and common name")

    # Hydrogenase reactions are missing protons
    write_to_log(output_log, f"Adding protons to hydrogenase reactions:")

    for rxn in tmodel.reactions:
        if 'HYD' in rxn.id:
            if 'ATP' not in rxn.id:
                write_to_log(output_log, f"- before: {rxn}")
                rxn.add_metabolites({tmodel.metabolites.h_e:2,
                                    tmodel.metabolites.h_c:-2})
                write_to_log(output_log, f"- after: {rxn}")

    # Remove duplicate/additional reactions
    write_to_log(output_log, f"Removing supplementary transport reactions: ")

    RXNS_TO_REMOVE = []

    for rxn in tmodel.reactions:
        if rxn.id.endswith("_add"):
            RXNS_TO_REMOVE.append(rxn)   
    tmodel.remove_reactions(RXNS_TO_REMOVE)

    write_to_log(output_log, f"Removed {RXNS_TO_REMOVE}")

    # Fix select reactions
    write_to_log(output_log, "Fixing select reactions FDH2, FDH3, TMAOR1e, TMAOR2e, GLCDe, DMSOR1e, DMSOR2e, SHCHF, SHCHD2 ")

    # FDH reactions are missing protons and water is on wrong side of membrane 
    tmodel.reactions.FDH2.add_metabolites({tmodel.metabolites.h_e:2,
                                        tmodel.metabolites.h_c:-2,
                                        tmodel.metabolites.h2o_c:1,
                                        tmodel.metabolites.h2o_e:-1})
    tmodel.reactions.FDH3.add_metabolites({tmodel.metabolites.h_e:2,
                                        tmodel.metabolites.h_c:-2,
                                        tmodel.metabolites.h2o_c:1,
                                        tmodel.metabolites.h2o_e:-1})

    # Missing transporterd metabolites
    tmodel.reactions.TMAOR1e.transported_mets = {tmodel.metabolites.tmao_e: -1}
    tmodel.reactions.TMAOR2e.transported_mets = {tmodel.metabolites.tmao_e: -1}
    tmodel.reactions.GLCDe.transported_mets = {tmodel.metabolites.get_by_id('glc-D_e'): -1}

    #D MSOR1e protons are incorrect - this resets the reaction to the bigg version 
    tmodel.reactions.DMSOR1e.add_metabolites({tmodel.metabolites.h_e:-2,
                                            tmodel.metabolites.h_c:2 })
    tmodel.reactions.DMSOR2e.add_metabolites({tmodel.metabolites.h_e:-2,
                                            tmodel.metabolites.h_c:2 })

    tmodel.reactions.DMSOR1e.transported_mets = {tmodel.metabolites.dmso_e: -1}
    tmodel.reactions.DMSOR2e.transported_mets = {tmodel.metabolites.dmso_e: -1}

    tmodel.reactions.SHCHF.add_metabolites({tmodel.metabolites.scl_c:-1,
                                            tmodel.metabolites.srch_c:1})

    # Sirohydrochlorin dehydrogenase is incorrectly defined in model 
    tmodel.reactions.SHCHD2.add_metabolites({tmodel.metabolites.scl_c:2,
                                            tmodel.metabolites.srch_c:-2})

    #Update srch metabolite name
    tmodel.metabolites.srch_c.annotation = {'bigg.metabolite': 'dscl'}

    # Specific reaction that invovled chemcial transformation as part of transport
    write_to_log(output_log, f"Defining transported metabolites of PTS reactions:")

    tmodel.reactions.NMNt7.transported_mets = {tmodel.metabolites.nmn_e:-1}

    # PTS mechanism 
    for rxn in tmodel.reactions:
        if 'pts' in rxn.id:
            rxn.transported_mets = {met:stoich for met, stoich in rxn.metabolites.items() if met.compartment == 'e'}
            write_to_log(output_log, f"{rxn}, {rxn.transported_mets}")

    #in imported model any charge metabolite represents a free cation that is transported 
    write_to_log(output_log, f"Defining transported charge of each reaction in the model")
    for rxn in tmodel.reactions:
        if 'biomass' not in rxn.id:   # ignore biomass reaction 
            
            transported_charge = {}
            for met, stoich in rxn.metabolites.items():
                if (met in tmodel.charge_dict.values()):
                    transported_charge[met.compartment] = stoich
            
            if len(transported_charge) != 0:
                rxn.transported_charge = transported_charge

    # Load default concentration bounds from the GAMS model:
    df_conc = hl.excel_to_df(gams)["ConcLimits"]
    write_to_log(output_log, f"Reading GAMS model file: {gams} (variable: 'ConcLimits')")

    # Rearrange data for easier use:
    df_conc = df_conc.reset_index()
    df_conc["met"] = df_conc["dim1"] + "_"+ df_conc["dim2"]
    df_conc = df_conc.pivot_table(columns="dim3", values="Value", index="met")

    write_to_log(output_log, f" - shape: {df_conc.shape}")
    write_to_log(output_log, f" - first 5 rows: {df_conc.head()}")

    # Change upper and lower concentration bounds the values in the dataframe df_conc (taken from GAMS)
    # (values are in mM)
    write_to_log(output_log, f"Setting concentration bounds:")
    for met, row in df_conc.iterrows():
        tmodel.metabolites.get_by_id(met).upper_bound = Q_(row["up"], "mM")
        tmodel.metabolites.get_by_id(met).lower_bound = Q_(row["lo"], "mM")
        write_to_log(output_log, f" - {met} ({tmodel.metabolites.get_by_id(met).lower_bound}, {tmodel.metabolites.get_by_id(met).upper_bound})")

    # Find correct biomass reaction
    biomass_rxns = [r for r in tmodel.reactions if 'biomass' in r.id.lower()]
    for r in biomass_rxns:
        print(r.id, ":", r.reaction[:80])

    # Correct reaction is 'biomass'

    write_to_log(output_log, f"Updating biomass metabolite (protons: H74; dfG0: -2.692234848 kJ/gDW) and adding growth-associated ATP maintenance:")

    # Define the protons in biomass 
    tmodel.metabolites.biomass_c.formula = 'H74'  # estimated http://www.ncbi.nlm.nih.gov/pubmed/20506321
    tmodel.metabolites.biomass_e.formula = 'H74'

    # Assign biomass
    tmodel.metabolites.biomass_c.biomass = True
    tmodel.metabolites.biomass_e.biomass = True

    tmodel.reactions.biomass.add_metabolites({tmodel.metabolites.atp_c: -31.2622,
                                                tmodel.metabolites.h2o_c: -31.2622,
                                                tmodel.metabolites.adp_c: +31.2622,
                                                tmodel.metabolites.pi_c:  +31.2622})


    # Define biomass formation energy: dfG0(biomass) [kJ gCDW-1] = -2.692234848 fom Battley 1991
    tmodel.metabolites.biomass_c.dfG0 = Q_(-2.692234848, "kJ/mol") * 1000    # values in J/gDW 
    tmodel.metabolites.biomass_e.dfG0 = Q_(-2.692234848, "kJ/mol") * 1000 

    tmodel.reactions.biomass_ce.ignore_snd = True
    write_to_log(output_log, f"Ignored 2nd law for biomass_ce transport reaction")

    write_to_log(output_log, f"Performing charge balance")
    for rxn in tmodel.reactions:
        thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)

    for met in tmodel.metabolites:
        if met.id in ['charge_c', 'charge_m', 'charge_e']:
            met.ignore_conc = True

    # Update thermodynamic information:
    write_to_log(output_log, f"Updating thermodynamic information")
    tmodel.update_thermo_info(fit_unknown_dfG0=True)

    return tmodel
