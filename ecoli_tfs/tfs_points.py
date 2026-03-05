import pandas as pd
import pta
import cobra
import cobra.io
import sys
import numpy as np
import glob

import tfs_scripts.preprocess
import tfs_scripts.main

IN_DIR = sys.argv[1]
IN_POINTS = sys.argv[2]
OUT_DIR = sys.argv[3]
NUM_SAMPLES_IN = int(sys.argv[4])
NUM_SAMPLES_OUT = int(sys.argv[5])

def create_tfs_model(tfs_dir, linear_constraints=False):

    tfs_files_dir = glob.glob(f"{tfs_dir}/*.*")

    model_file = [x for x in tfs_files_dir if x.endswith(".sbml")][0].strip()
    #vbounds = [x for x in tfs_files_dir if x.endswith("vbounds.csv")][0].strip() # Bounds are already applied to loaded model
    lcv = [x for x in tfs_files_dir if x.endswith("lcv.csv")][0].strip()
    lcm = [x for x in tfs_files_dir if x.endswith("lcm.csv")][0].strip()
    lcb = [x for x in tfs_files_dir if x.endswith("lcb.csv")][0].strip()
    drg0pm = [x for x in tfs_files_dir if x.endswith("drg0pm.csv")][0].strip()
    drg0cs = [x for x in tfs_files_dir if x.endswith("drg0cs.csv")][0].strip()
    restrained_rxns_csv = [x for x in tfs_files_dir if x.endswith("ignored_rxns.csv")][0].strip()

    restrained = pd.read_csv(restrained_rxns_csv, index_col=0)
    restrained_rxns = []
    for x in restrained.itertuples():
        restrained_rxns.append(x.rxn)

    print(restrained_rxns)

    tfs_model = tfs_scripts.main.main(model_file, None, lcm, lcv, lcb, drg0pm, drg0cs, linear_constraitns=linear_constraints, restrained_rxns=restrained_rxns)

    return tfs_model
    
def save_drg_results_for_cluster(tfs_model, tfs_result, outfolder):
    #cobra.io.write_sbml_model(tfs_model.model, f"{outfolder}/model.sbml")
    tfs_result.samples.to_csv(f"{outfolder}/samples.csv")
    tfs_result.orthants.to_csv(f"{outfolder}/orthants.csv")
    tfs_result.psrf.to_csv(f"{outfolder}/psrf.csv")

#points = pd.read_csv(IN_POINTS,index_col=0).to_numpy()

tfs_model = create_tfs_model(IN_DIR, True)

points = 0
while points is not None:
  try:
    points = tfs_model.get_initial_points(NUM_SAMPLES_IN)
  except Exception as e:
    print("Skipped infeasible point")
  
  if points is not None:
    df = pd.DataFrame(points)
    df.to_csv(f"{OUT_DIR}/points_{NUM_SAMPLES_IN}.csv")
    break

#result = pta.sample_drg(tfs_model, initial_points=points, max_threads=16)
#save_drg_results_for_cluster(tfs_model, result, OUT_DIR)

#drG_samples = result.samples.loc[np.random.choice(result.samples.index, NUM_SAMPLES_IN, replace=False)]
#fluxes = pta.sample_fluxes_from_drg(tfs_model.model, drG_samples, result.orthants, num_approx_samples=NUM_SAMPLES_OUT)
#fluxes.to_feather(f"{OUT_DIR}/fluxes_{NUM_SAMPLES_IN}_{NUM_SAMPLES_OUT}.feather")