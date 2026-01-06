
import os.path as path
import glob, os
from datetime import datetime
from importlib.metadata import version
import cobra
import thermo_flux
from thermo_flux.io import load_excel as ex
from thermo_flux.core.model import ThermoModel
from thermo_flux.tools.drg_tools import calc_dfG_transform
from equilibrator_api import  Q_
import pandas as pd
from thermo_flux.io import helper_load as hl
import numpy as np
from thermo_flux.io import load_excel as ex
import gurobipy as gp
from scripts.logger import write_to_log
from gurobipy import GRB
from cobra.flux_analysis import flux_variability_analysis

CONDITIONS_REGRESS = ["WT-Glc_I", "WT-Gal_I", "WT-Fruc_I", "WT-Mann_I", "dptsG-Glc_I", "WT-Ace_I", "WT-Succ_I", "WT-Fum_I", "WT-Glyc_I", "WT-Pyr_I", "WT-GlyCAA_II"];

def gen_model_from_core(name: str, model_data: str, kegg: str, reed: str, inchi: str, gams: str, add_o2: bool, add_co2: bool, update_thermodynamics = True):
    model = cobra.io.load_model('e_coli_core')
    tmodel = ThermoModel(model)

    # # ADD HYDROXYBENZOATE TRANSPORT & EXCHANGE:

    # # Define extracellular orotate:
    # hbz_e = cobra.Metabolite(id="4hbz_e", compartment="e")
    # hbz_e = thermo_flux.core.metabolite.ThermoMetabolite(hbz_e, model=tmodel)
    # hbz_e.annotation = tmodel.metabolites.get_by_id("4hbz_c").annotation  
    # tmodel.metabolites.append(hbz_e)

    # # Hydroxybenzoate H+ antiporter (aaeB):
    # HBZt3 = cobra.Reaction("HBZt3")
    # HBZt3.lower_bound = -1000
    # HBZt3.upper_bound = 1000
    # HBZt3.add_metabolites({tmodel.metabolites.get_by_id("4hbz_e"): +1,
    #                   tmodel.metabolites.h_e: -1,
    #                   tmodel.metabolites.get_by_id("4hbz_c"): -1,
    #                   tmodel.metabolites.h_c: +1})
    # HBZt3 = thermo_flux.core.reaction.ThermoReaction(HBZt3, model=tmodel)

    # # Exchange:
    # EX_4hbz = cobra.Reaction("EX_4hbz")
    # EX_4hbz.add_metabolites({tmodel.metabolites.get_by_id("4hbz_e"): -1})
    # EX_4hbz = thermo_flux.core.reaction.ThermoReaction(EX_4hbz, model=tmodel)

    # # Add reactions to the model
    # tmodel.add_reactions([HBZt3, EX_4hbz])

    # # ADD ISOPROPYLMALATE TRANSPORT & EXCHANGE:

    # # Define extracellular orotate:
    # ipm_e = cobra.Metabolite(id="3c3hmp_e", compartment="e")
    # ipm_e = thermo_flux.core.metabolite.ThermoMetabolite(ipm_e, model=tmodel)
    # ipm_e.annotation = tmodel.metabolites.get_by_id("3c3hmp_c").annotation  
    # tmodel.metabolites.append(ipm_e)

    # # Diffusion across the membrane:
    # ipm_diff = cobra.Reaction("IPMex")
    # ipm_diff.lower_bound = -1000
    # ipm_diff.upper_bound = 1000
    # ipm_diff.add_metabolites({tmodel.metabolites.get_by_id("3c3hmp_c"): -1,
    #                         tmodel.metabolites.get_by_id("3c3hmp_e"): +1})
    # ipm_diff = thermo_flux.core.reaction.ThermoReaction(ipm_diff, model=tmodel)


    # # Exchange:
    # EX_3c3hmp = cobra.Reaction("EX_3c3hmp")
    # EX_3c3hmp.add_metabolites({tmodel.metabolites.get_by_id("3c3hmp_e"): -1})
    # EX_3c3hmp = thermo_flux.core.reaction.ThermoReaction(EX_3c3hmp, model=tmodel)

    # # Add reactions to the model
    # tmodel.add_reactions([ipm_diff, EX_3c3hmp])


    # # Define extracellular orotate:
    # oro_e = cobra.Metabolite(id="orot_e", compartment="e")
    # oro_e = thermo_flux.core.metabolite.ThermoMetabolite(oro_e, model=tmodel)
    # oro_e.annotation = tmodel.metabolites.orot_c.annotation  
    # tmodel.metabolites.append(oro_e)

    # # Diffusion across the membrane:
    # oro_diff = cobra.Reaction("OROTex")
    # oro_diff.lower_bound = -1000
    # oro_diff.upper_bound = 1000
    # oro_diff.add_metabolites({tmodel.metabolites.orot_c: -1,
    #                         tmodel.metabolites.orot_e: +1})
    # oro_diff = thermo_flux.core.reaction.ThermoReaction(oro_diff, model=tmodel)

    # # Dicarboxylate/H+ symporter (dctA):
    # dcta = cobra.Reaction("DCTA")
    # dcta.lower_bound = -1000
    # dcta.upper_bound = 0
    # dcta.add_metabolites({tmodel.metabolites.orot_e: +1,
    #                     tmodel.metabolites.h_e: +1,
    #                     tmodel.metabolites.orot_c: -1,
    #                     tmodel.metabolites.h_c: -1})
    # dcta = thermo_flux.core.reaction.ThermoReaction(dcta, model=tmodel)

    # # Exchange:
    # EX_oro = cobra.Reaction("EX_oro")
    # EX_oro.add_metabolites({tmodel.metabolites.orot_e: -1})
    # EX_oro = thermo_flux.core.reaction.ThermoReaction(EX_oro, model=tmodel)

    # # Add reactions to the model
    # tmodel.add_reactions([oro_diff, dcta, EX_oro])

    tmodel.pH = {"c": Q_(7.6), "e": Q_(7)} #pH
    tmodel.I = {"c": Q_(0.25,'M'), "e": Q_(0.25,'M')} #ionic stength
    tmodel.phi = {'ce':Q_(0.15,'V')} #membrane potential ‘ce’ represents the voltage between compartment ‘c’ and compartment 'e’ defined as Phic - Phie
    tmodel.pMg = {'e': Q_(3), 'c': Q_(3)}

    for met in tmodel.metabolites:
        met.annotation["bigg.metabolite"] = met.id[:-2]


    # Update the inchi strings of some unknown metabolites
    df = pd.read_csv(reed, header=None).set_index(0)

    unknown_mets = []
    for met in tmodel.metabolites:
        if met.id[:-2] in df.index:
            met.annotation['InChI'] = df.loc[met.id[:-2]].iloc[1]
            unknown_mets.append(met)

    for met in tmodel.metabolites:
        if ('kegg' in met.annotation):
            if (met.annotation['kegg'] in df.index):
            
                inchi = df.loc[met.annotation['kegg']].iloc[0]

                if type(inchi) is str:
                    met.annotation['InChI'] = df.loc[met.annotation['kegg']].iloc[0]
                    unknown_mets.append(met)

    # Additional data from excel spredsheet needs to be imported for some compounds with unknown structure
    sheets = ex.read_excelfile(model_data)

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

    
    for rxn in tmodel.reactions:
        if 'HYD' in rxn.id:
            if 'ATP' not in rxn.id:
                rxn.add_metabolites({tmodel.metabolites.h_e:2,
                                    tmodel.metabolites.h_c:-2})


    RXNS_TO_REMOVE = []

    for rxn in tmodel.reactions:
        if rxn.id.endswith("_add"):
            RXNS_TO_REMOVE.append(rxn)
            print(f"Removing: {rxn.id}")
    tmodel.remove_reactions(RXNS_TO_REMOVE)

    # FDH reactions are missing protons and water is on wrong side of membrane 
    # tmodel.reactions.FDH2.add_metabolites({tmodel.metabolites.h_e:2,
    #                                     tmodel.metabolites.h_c:-2,
    #                                     tmodel.metabolites.h2o_c:1,
    #                                     tmodel.metabolites.h2o_e:-1})
    # tmodel.reactions.FDH3.add_metabolites({tmodel.metabolites.h_e:2,
    #                                     tmodel.metabolites.h_c:-2,
    #                                     tmodel.metabolites.h2o_c:1,
    #                                     tmodel.metabolites.h2o_e:-1})

    # # Missing transporterd metabolites
    # tmodel.reactions.TMAOR1e.transported_mets = {tmodel.metabolites.tmao_e: -1}
    # tmodel.reactions.TMAOR2e.transported_mets = {tmodel.metabolites.tmao_e: -1}
    # tmodel.reactions.GLCDe.transported_mets = {tmodel.metabolites.get_by_id('glc-D_e'): -1}

    # #D MSOR1e protons are incorrect - this resets the reaction to the bigg version 
    # tmodel.reactions.DMSOR1e.add_metabolites({tmodel.metabolites.h_e:-2,
    #                                         tmodel.metabolites.h_c:2 })
    # tmodel.reactions.DMSOR2e.add_metabolites({tmodel.metabolites.h_e:-2,
    #                                         tmodel.metabolites.h_c:2 })

    # tmodel.reactions.DMSOR1e.transported_mets = {tmodel.metabolites.dmso_e: -1}
    # tmodel.reactions.DMSOR2e.transported_mets = {tmodel.metabolites.dmso_e: -1}

    # tmodel.reactions.SHCHF.add_metabolites({tmodel.metabolites.scl_c:-1,
    #                                         tmodel.metabolites.srch_c:1})

    # # Sirohydrochlorin dehydrogenase is incorrectly defined in model 
    # tmodel.reactions.SHCHD2.add_metabolites({tmodel.metabolites.scl_c:2,
    #                                         tmodel.metabolites.srch_c:-2})

    # #Update srch metabolite name
    # tmodel.metabolites.srch_c.annotation = {'bigg.metabolite': 'dscl'}

    # tmodel.reactions.NMNt7.transported_mets = {tmodel.metabolites.nmn_e:-1}

    # PTS mechanism 
    for rxn in tmodel.reactions:
        if 'pts' in rxn.id:
            rxn.transported_mets = {met:stoich for met, stoich in rxn.metabolites.items() if met.compartment == 'e'}

    #in imported model any charge metabolite represents a free cation that is transported 
    #for rxn in tmodel.reactions:
    #    if 'biomass' not in rxn.id:   # ignore biomass reaction 
    #        
    #        transported_charge = {}
    #        for met, stoich in rxn.metabolites.items():
    #            if (met in tmodel.charge_dict.values()):
    #                transported_charge[met.compartment] = stoich
    #        
    #        if len(transported_charge) != 0:
    #            rxn.transported_charge = transported_charge

    # Load default concentration bounds from the GAMS model:
    #df_conc = hl.excel_to_df(gams)["ConcLimits"]

    # Rearrange data for easier use:
    #df_conc = df_conc.reset_index()
    #df_conc["met"] = df_conc["dim1"] + "_"+ df_conc["dim2"]
    #df_conc = df_conc.pivot_table(columns="dim3", values="Value", index="met")

    # Change upper and lower concentration bounds the values in the dataframe df_conc (taken from GAMS)
    # (values are in mM)
    #met_ids = [m.id for m in tmodel.metabolites]

    #for met, row in df_conc.iterrows():
    #    if met in met_ids:
    #        tmodel.metabolites.get_by_id(met).upper_bound = Q_(row["up"], "mM")
    #        tmodel.metabolites.get_by_id(met).lower_bound = Q_(row["lo"], "mM")

    # Find correct biomass reaction
    biomass_rxns = [r for r in tmodel.reactions if 'biomass' in r.id.lower()]
    for r in biomass_rxns:
        print(r.id, ":", r.reaction[:80])

    # Correct reaction is 'biomass'
    # Define the protons in biomass 
    #tmodel.metabolites.biomass_c.formula = 'H74'  # estimated http://www.ncbi.nlm.nih.gov/pubmed/20506321
    #tmodel.metabolites.biomass_e.formula = 'H74'

    # Assign biomass
    #tmodel.metabolites.biomass_c.biomass = True
    #tmodel.metabolites.biomass_e.biomass = True

    #tmodel.reactions.BIOMASS_Ecoli_core_w_GAM.add_metabolites({tmodel.metabolites.atp_c: -31.2622,
    #                                            tmodel.metabolites.h2o_c: -31.2622,
    #                                            tmodel.metabolites.adp_c: +31.2622,
    #                                            tmodel.metabolites.pi_c:  +31.2622})


    # Define biomass formation energy: dfG0(biomass) [kJ gCDW-1] = -2.692234848 fom Battley 1991
    # structural UCFW = 24.190 Da (Battley 1991)
    # dGf biomass: −65.10 kJ/C-mol (Battley 1991, referenced in DOI:10.3390/e22030277)
    # dGf per gCDW = -2.69119470856

    #energy_diff = calc_dfG_transform(tmodel.metabolites.biomass_e) - calc_dfG_transform(tmodel.metabolites.biomass_c)

    #base_dfg = (Q_(-2.692234848, "kJ/mol") * 1000) # values in J/gDW 

    #tmodel.metabolites.biomass_c.dfG0 = base_dfg
    #tmodel.metabolites.biomass_e.dfG0 = base_dfg

    #tmodel.reactions.biomass_ce.ignore_snd = True

    if update_thermodynamics:
        for rxn in tmodel.reactions:
            thermo_flux.tools.drg_tools.reaction_balance(rxn, balance_charge=True, balance_mg=False)

        for met in tmodel.metabolites:
            if met.id in ['charge_c', 'charge_m', 'charge_e']:
                met.ignore_conc = True

        # Update thermodynamic information:
        tmodel.update_thermo_info(fit_unknown_dfG0=True)

    return tmodel

    

def gen_model(name: str, model_xlsx: str, kegg: str, reed: str, inchi:str, gams: str, output_log: str, add_o2: bool, add_co2: bool, update_thermodynamics=True):
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

    # ADD OROTATE TRANSPORT AND EXCHANGE

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
    # structural UCFW = 24.190 Da (Battley 1991)
    # dGf biomass: −65.10 kJ/C-mol (Battley 1991, referenced in DOI:10.3390/e22030277)
    # dGf per gCDW = -2.69119470856

    #energy_diff = calc_dfG_transform(tmodel.metabolites.biomass_e) - calc_dfG_transform(tmodel.metabolites.biomass_c)

    base_dfg = (Q_(-2.692234848, "kJ/mol") * 1000) # values in J/gDW 

    tmodel.metabolites.biomass_c.dfG0 = base_dfg
    tmodel.metabolites.biomass_e.dfG0 = base_dfg

    tmodel.reactions.biomass_ce.ignore_snd = True
    write_to_log(output_log, f"Ignored 2nd law for biomass_ce transport reaction")

    if update_thermodynamics:
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


def apply_physio_data(tmodel, condition :str, input_exp: str, input_conc: str, input_metabolomics: str, input_gams: str, relax_flux_bounds, include_CO2: bool, include_O2: bool, allow_other_excr: bool, output_log: str, open_exchanges=False):
    "Apply regression to the Tmodel (not gurobi model) for FBA and blocked reaction analysis"
    df_conc = hl.excel_to_df(input_gams)["ConcLimits"]

    # Rearrange data for easier use:
    df_conc = df_conc.reset_index()
    df_conc["met"] = df_conc["dim1"] + "_"+ df_conc["dim2"]
    df_conc = df_conc.pivot_table(columns="dim3", values="Value", index="met")

    # Import experimental data:
    reg_data = pd.read_csv(input_exp)
    write_to_log(output_log, f"Reading experimental flux data: {input_exp}")

    reg_data.set_index(["cond", "rxn"], inplace=True) 
    reg_data.head()

    # Store gas fluxes:
    reg_data_gas = reg_data.swaplevel().copy()
    reg_data_gas = reg_data_gas.loc[["EX_co2", "EX_o2"]]
    reg_data_gas = reg_data_gas.swaplevel()
    reg_data_gas

    if include_CO2 is False:
        reg_data_no_gas = reg_data.swaplevel().copy()
        reg_data_no_gas = reg_data_no_gas.drop(["EX_co2"])
        reg_data_no_gas = reg_data_no_gas.swaplevel()
        reg_data = reg_data_no_gas
        write_to_log(output_log, f" - ignoring CO2 data")
        
        
    if include_O2 is False:
        reg_data_no_gas = reg_data.swaplevel().copy()
        reg_data_no_gas = reg_data_no_gas.drop(["EX_o2"]) 
        reg_data_no_gas = reg_data_no_gas.swaplevel()
        reg_data = reg_data_no_gas
        write_to_log(output_log, f" - ignoring O2 data")

    # Import experimental data:
    conc_data = pd.read_csv(input_conc)
    write_to_log(output_log, f"Reading experimental extracellular concentration data: {input_conc}")

    conc_data.set_index(["cond", "met"], inplace=True) 
    conc_data.head()

    #Should probably change this to just the given condition..
    volume_data = pd.DataFrame({cond: {"c": 1.0} for cond in CONDITIONS_REGRESS} ).T #specify the volume fractions for each condition
    volume_data.head()

    # Import experimental data:
    write_to_log(output_log, f"Reading intracellular metabolite concentration data: {input_metabolomics}")
    met_data = pd.read_csv(input_metabolomics)
    met_data.set_index(["cond", "met"], inplace=True) 
    met_data.head()

    conds_with_data = list(met_data.reset_index().cond.unique())
    missing_conds = [cond for cond in CONDITIONS_REGRESS if cond not in conds_with_data]

    df_missing = pd.DataFrame({"cond": missing_conds, 
                            "met": "g6p",   # code doesn't deal well if we type a non-existince met here...
                            "mean": np.nan, 
                            "sd": np.nan, }).set_index(["cond", "met"])

    met_data = pd.concat([met_data, df_missing])

    # Store the indices of all reactions:
    map_rxn_id = {rxn.id: index for index, rxn in enumerate(tmodel.reactions)}

    exchanges = [rxn.id for rxn in tmodel.exchanges]
    print(exchanges)

    exchanges_to_relax = ["EX_C", "EX_h", "EX_h2o", "EX_k", "EX_nh3", "EX_pi", "EX_so4"]

    if include_CO2 is False:
        exchanges_to_relax += ["EX_co2"]
        
    if include_O2 is False:
        exchanges_to_relax += ["EX_o2"]

    if allow_other_excr is True:
        upper_bound_exchanges = 100
    else:
        upper_bound_exchanges = 0

    settings_tfba = {"error_type": "covariance",
                    "qnorm": 1,
                    "alpha": 0.95, 
                    "epsilon": 0.5,
                    "nullspace": None,
                    "gdiss_constraint": True,
                    "sigmac_limit": 130}


    settings_regression = {"flux_data": reg_data,
                        "metabolite_data": met_data,
                        "volume_data": volume_data,
                        "conc_units": "mM",
                        "conc_fit": False,
                        "flux_fit": True,
                        "drG_fit": True, 
                        "resnorm": 1, 
                        "error_type": "covariance"}

    write_to_log(output_log, f"Setting up regressions:")
    write_to_log(output_log, f" - exchanges to be relaxed: {exchanges_to_relax}")
    write_to_log(output_log, f" - stdev of experimental fluxes increased by factor of: {relax_flux_bounds}")
    write_to_log(output_log, f" - settings for tFBA: {settings_tfba}")
    write_to_log(output_log, f" - settings for regression: {settings_regression}")

    # Reset all flux bounds to +- 100:
    for rxn in tmodel.reactions:
        tmodel.reactions.get_by_id(rxn.id).lower_bound = -100
        tmodel.reactions.get_by_id(rxn.id).upper_bound = 100
    
    # Add non-growth associate ATP maintenance cost:
    tmodel.reactions.ATPHYD.lower_bound = 3.15

    # Fix exchange reaction directions:
    for rxn in exchanges:
        tmodel.reactions.get_by_id(rxn).lower_bound = 0
        tmodel.reactions.get_by_id(rxn).upper_bound = upper_bound_exchanges

    # Allow for excretion of selected metabolites:
    #for rxn in EXCEPTIONS:
    #    tmodel.reactions.get_by_id(rxn).lower_bound = 0
    #    tmodel.reactions.get_by_id(rxn).upper_bound = +100
        
    # Relax essential exchanges:
    for rxn_rel in exchanges_to_relax:
        tmodel.reactions.get_by_id(rxn_rel).lower_bound = -100
        tmodel.reactions.get_by_id(rxn_rel).upper_bound = +100

    # Only allow secretion for carbon sources other than the current condition
    # Gets all exchange reactions for every condition and sets it to export only
    # Export reactions for this specific condition will still be set to the physiological data in the step below
    #if other_carbon_export_only:   
        #exp_data = pd.read_csv(input_exp) 
        #carbon_source_exchanges = [row["rxn"] for _, row in exp_data.iterrows()]
        #print(carbon_source_exchanges)
        
        #CONDITIONS_TO_REGRESS = ["WT-Glc_I", "WT-Gal_I", "WT-Fruc_I", "WT-Mann_I", "dptsG-Glc_I", 
#                         "WT-Ace_I", "WT-Succ_I", "WT-Fum_I", "WT-Glyc_I", "WT-Pyr_I",
#                         "WT-GlyCAA_II"]

     #   carbon_source_exchanges =["EX_glc", "EX_gal", "EX_fru", "EX_man", "EX_ac", "EX_succ", "EX_fum", "EX_glyc", "EX_pyr", "EX_gly"]
        
     #   exchange_rxn_set = set(carbon_source_exchanges)
     #   for rxn in exchange_rxn_set:
     #       tmodel.reactions.get_by_id(rxn).lower_bound = 0
     #       tmodel.reactions.get_by_id(rxn).upper_bound = 100

    # Fix flux for the measured exchange reactions:
    for rxn, row in reg_data.loc[condition].iterrows():
        tmodel.reactions.get_by_id(rxn).lower_bound = -100
        tmodel.reactions.get_by_id(rxn).upper_bound = 100

        if not open_exchanges:
            tmodel.reactions.get_by_id(rxn).lower_bound = row["mean"] - relax_flux_bounds * row["sd"]
            tmodel.reactions.get_by_id(rxn).upper_bound = row["mean"] + relax_flux_bounds * row["sd"]
            #write_to_log(output_log, f" - {rxn}: ({tmodel.reactions.get_by_id(rxn).lower_bound :.3}, {tmodel.reactions.get_by_id(rxn).upper_bound :.3})")

    if condition.startswith("dptsG-Glc"):
        tmodel.reactions.GLCpts.lower_bound = 0
        tmodel.reactions.GLCpts.upper_bound = 0
        write_to_log(output_log, f" - blocked GLCpts")
        
    # Set metabolite concentrations to the values in the GAMS model:   
    for met, row in df_conc.iterrows():
        tmodel.metabolites.get_by_id(met).upper_bound = Q_(row["up"], "mM")
        tmodel.metabolites.get_by_id(met).lower_bound = Q_(row["lo"], "mM")
        
    # Fix concentration for the measured extracellular metabolites:
    for met, row in conc_data.loc[condition].iterrows():
        tmodel.metabolites.get_by_id(met).lower_bound = Q_(1e-9, "M")
        tmodel.metabolites.get_by_id(met).upper_bound = Q_(100, "M")
        tmodel.metabolites.get_by_id(met).lower_bound = Q_(row["conc_M_min"], "M")
        tmodel.metabolites.get_by_id(met).upper_bound = Q_(row["conc_M_max"], "M")
        write_to_log(output_log, f" - {met}: ({tmodel.metabolites.get_by_id(met).lower_bound :.3}, {tmodel.metabolites.get_by_id(met).upper_bound :.3})")

    return tmodel

def constrain_bounds_fva(tmodel, output_log: str):
    "Tighten bounds of a ThermoModel's reactions"
    all_reactions = [r.id for r in tmodel.reactions]
    fva = flux_variability_analysis(tmodel, reaction_list=all_reactions, fraction_of_optimum=1.0, processes=1)

    for r_id, row in fva.iterrows():
        if r_id not in tmodel.reactions:
            write_to_log(output_log, f"{r_id} is not a valid reaction. Skipping.")
            continue
        r = tmodel.reactions.get_by_id(r_id)
        fmin, fmax = float(row["minimum"]), float(row["maximum"])

        write_to_log(output_log, f"FVA: {r.id} bounds are {fmin}, {fmax}. Experimental bounds are {r.lower_bound}, {r.upper_bound} ")

        lower_bound = max(r.lower_bound, fmin)
        upper_bound = min(r.upper_bound, fmax)

        if lower_bound > upper_bound:
            lower_bound = fmin
            upper_bound = fmax
            write_to_log(output_log, f"Lower bound higher than upper ({lower_bound}, {upper_bound}). Falling back to {fmin}, {fmax}")
        
        write_to_log(output_log, f"{r.id} new bounds are {lower_bound}, {upper_bound}")

        r.lower_bound, r.upper_bound = lower_bound, upper_bound

def run_fva(tmodel, output_log: str):
    "Tighten bounds of a ThermoModel's reactions"
    bounds = []
    all_reactions = [r.id for r in tmodel.reactions]
    fva = flux_variability_analysis(tmodel, reaction_list=all_reactions, fraction_of_optimum=1.0, processes=1)

    for r_id, row in fva.iterrows():
        if r_id not in tmodel.reactions:
            write_to_log(output_log, f"{r_id} is not a valid reaction. Skipping.")
            continue
        r = tmodel.reactions.get_by_id(r_id)
        fmin, fmax = float(row["minimum"]), float(row["maximum"])

        write_to_log(output_log, f"FVA: {r.id} bounds are {fmin}, {fmax}. Experimental bounds are {r.lower_bound}, {r.upper_bound} ")

        lower_bound = max(r.lower_bound, fmin)
        upper_bound = min(r.upper_bound, fmax)

        if lower_bound > upper_bound:
            lower_bound = fmin
            upper_bound = fmax
            write_to_log(output_log, f"Lower bound higher than upper ({lower_bound}, {upper_bound}). Falling back to {fmin}, {fmax}")
        
        write_to_log(output_log, f"{r.id} new bounds are {lower_bound}, {upper_bound}")

        #r.lower_bound, r.upper_bound = lower_bound, upper_bound
        bounds.append((lower_bound, upper_bound))
    return bounds

def run_optimization(tmodel, name: str, condition :str, input_exp: str, input_conc: str, input_metabolomics: str, input_gams: str, output_log: str, include_CO2: bool, include_O2: bool, allow_other_excr: bool):

    df_conc = hl.excel_to_df(input_gams)["ConcLimits"]

    # Rearrange data for easier use:
    df_conc = df_conc.reset_index()
    df_conc["met"] = df_conc["dim1"] + "_"+ df_conc["dim2"]
    df_conc = df_conc.pivot_table(columns="dim3", values="Value", index="met")

    # Import experimental data:
    reg_data = pd.read_csv(input_exp)
    write_to_log(output_log, f"Reading experimental flux data: {input_exp}")

    reg_data.set_index(["cond", "rxn"], inplace=True) 
    #reg_data.head()

    # Store gas fluxes:
    reg_data_gas = reg_data.swaplevel().copy()
    reg_data_gas = reg_data_gas.loc[["EX_co2", "EX_o2"]]
    reg_data_gas = reg_data_gas.swaplevel()
    #reg_data_gas

    if include_CO2 is False:
        reg_data_no_gas = reg_data.swaplevel().copy()
        reg_data_no_gas = reg_data_no_gas.drop(["EX_co2"])
        reg_data_no_gas = reg_data_no_gas.swaplevel()
        reg_data = reg_data_no_gas
        write_to_log(output_log, f" - ignoring CO2 data")
        
        
    if include_O2 is False:
        reg_data_no_gas = reg_data.swaplevel().copy()
        reg_data_no_gas = reg_data_no_gas.drop(["EX_o2"]) 
        reg_data_no_gas = reg_data_no_gas.swaplevel()
        reg_data = reg_data_no_gas
        write_to_log(output_log, f" - ignoring O2 data")

    # Import experimental data:
    conc_data = pd.read_csv(input_conc)
    write_to_log(output_log, f"Reading experimental extracellular concentration data: {input_conc}")

    conc_data.set_index(["cond", "met"], inplace=True) 
    conc_data.head()

    #Should probably change this to just the given condition..
    volume_data = pd.DataFrame({cond: {"c": 1.0} for cond in CONDITIONS_REGRESS} ).T #specify the volume fractions for each condition
    volume_data.head()

    # Import experimental data:
    write_to_log(output_log, f"Reading intracellular metabolite concentration data: {input_metabolomics}")
    met_data = pd.read_csv(input_metabolomics)
    met_data.set_index(["cond", "met"], inplace=True) 
    met_data.head()

    conds_with_data = list(met_data.reset_index().cond.unique())
    missing_conds = [cond for cond in CONDITIONS_REGRESS if cond not in conds_with_data]

    df_missing = pd.DataFrame({"cond": missing_conds, 
                            "met": "g6p",   # code doesn't deal well if we type a non-existince met here...
                            "mean": np.nan, 
                            "sd": np.nan, }).set_index(["cond", "met"])

    met_data = pd.concat([met_data, df_missing])

    # Store the indices of all reactions:
    map_rxn_id = {rxn.id: index for index, rxn in enumerate(tmodel.reactions)}

    exchanges = [rxn.id for rxn in tmodel.exchanges]

    exchanges_to_relax = ["EX_C", "EX_h", "EX_h2o", "EX_k", "EX_nh3", "EX_pi", "EX_so4"]

    if include_CO2 is False:
        exchanges_to_relax += ["EX_co2"]
        
    if include_O2 is False:
        exchanges_to_relax += ["EX_o2"]

    if allow_other_excr is True:
        upper_bound_exchanges = 100
    else:
        upper_bound_exchanges = 0

    settings_tfba = {"error_type": "covariance",
                    "qnorm": 1,
                    "alpha": 0.95, 
                    "epsilon": 0.5,
                    "nullspace": None,
                    "gdiss_constraint": True,
                    "sigmac_limit": 130}


    settings_regression = {"flux_data": reg_data,
                        "metabolite_data": met_data,
                        "volume_data": volume_data,
                        "conc_units": "mM",
                        "conc_fit": False,
                        "flux_fit": True,
                        "drG_fit": True, 
                        "resnorm": 1, 
                        "error_type": "covariance"}

    #Initialize model:
    tmodel.m = None  # clear any previously build optimization models 
    tmodel.objective = tmodel.reactions.biomass_EX      # needed, otherwise add_TFBA_variables() gives an error due to lack of an obj. function

    # Add thermodynamics constraints
    tmodel.add_TFBA_variables(conds=[condition], **settings_tfba) 

    # Setup regression to experimental data:
    tmodel.regression([condition], **settings_regression)
    tmodel.m.update()

    
    # Display objective prior to adding RQ as an additional objective:
#     print("Before: ", tmodel.m.getObjective())
    write_to_log(output_log, f" - Objective function before: {tmodel.m.getObjective()}")
    
    
    # Determine RQ for the given condition:
    data_gas = reg_data_gas.loc[condition]
    vo2, vo2_err = data_gas.loc["EX_o2"]
    vco2, vco2_err = data_gas.loc["EX_co2"]
    rq = - vco2 / vo2
    rq_err = rq * np.sqrt( (vo2_err / vo2)**2 + (vco2_err / vco2)**2)
    write_to_log(output_log, f" - RQ: {rq :.2} (vCO2 = {vco2 :.3} / vO2 = {vo2 :.3})")

    
    # Add new constraint - residual of RQ: 
    resrq = tmodel.m.addMVar(lb=0, ub=GRB.INFINITY, shape=(1,1), name="resRQ") 
    mrq = tmodel.m.addMVar(lb=0, ub=GRB.INFINITY, shape=(1,1), name="RQ") 
    tmodel.mvars["resRQ"] = resrq
    tmodel.mvars["rq"] = mrq
    
    tmodel.m.addConstr(tmodel.mvars["resRQ"][0, 0] >= ( (mrq - rq)/rq_err ), name = "resRQ_pos")
    tmodel.m.addConstr(tmodel.mvars["resRQ"][0, 0] >= (-(mrq - rq)/rq_err ), name = "resRQ_neg")
    
    
    # Impose RQ on the gas fluxes:
    idx_o2 = map_rxn_id["EX_o2"]
    idx_co2 = map_rxn_id["EX_co2"]   
    tmodel.m.addConstr(tmodel.mvars["v"][0, idx_co2] == (-mrq)*tmodel.mvars["v"][0, idx_o2], name="enforce_RQ")
    tmodel.m.update()

    # Update objective function to include the new RQ constraint:
    # (+ adjust weight of each term, giving equal weight to the RQ and each of the fluxes being fit)
    no_fluxes = len(list(reg_data.unstack()["mean"].T.index) )
    initial_weight = 1/no_fluxes
    new_weight = 1/(no_fluxes + 1)   # to account for RQ constraint
    
    # Equal weight to RQ and resflx:
    tmodel.m.setObjective(tmodel.m.getObjective() - initial_weight*tmodel.mvars["resflx"].sum() + new_weight*(resrq.sum() + tmodel.mvars["resflx"].sum()), GRB.MINIMIZE)
    
    write_to_log(output_log, f" - Objective function after: {tmodel.m.getObjective()}")
    
    # Write model for optimization on the HPC cluster:
    # model_filename = f"{output_log}{path.sep}single_{CONDITION}_fit-co2-{INCLUDE_CO2}_fit-o2-{INCLUDE_O2}_allow-excretion-{ALLOW_OTHER_EXCRETION}.mps"
    # tmodel.m.write(model_filename)
    # write_to_log(output_log, f"Saved model to: {model_filename}")

    tmodel.m.Params.TimeLimit = 3600
    tmodel.m.Params.Threads = 16    
    tmodel.m.optimize()
    
    try:
        sols = tmodel.solution()
        best_solution = sols[-1]
        best_solution.to_csv(f"solutions{path.sep}{name}_{condition}_SOLUTION.csv")
    except:
        print("Failed to save solution to csv, falling back to Gurobi .sol output")
        tmodel.m.write(f"solutions{path.sep}{name}_{condition}_SOLUTION.sol")

    return tmodel