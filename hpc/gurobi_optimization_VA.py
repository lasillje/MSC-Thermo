import sys
##modify path so we get gurobipy from venv to use gp.read
import os

sys.path.append('/home4/p313427/freeflux/') ## make cleaner with user expand
sys.path.append('/home4/p313427/')
sys.path.append(os.path.abspath('../'))
sys.path.append(os.path.abspath('../sampler/'))
sys.path.append('/home4/p313427/sampler/')
sys.path.append('/home4/p313427/thermo_flux/')

from gurobipy import read, GRB


print(sys.argv)

# Parse the input arguments:
FILE = str(sys.argv[1])
FNAME=sys.argv[2]
HOURS_OPTIM = int(sys.argv[3])
CPUS_OPTIM = int(sys.argv[4])
OUTPUT_FOLDER=sys.argv[5]

print(
    f"Filename: {FILE}\nDuration of optimization: {HOURS_OPTIM} h\nCPUs used: {CPUS_OPTIM}"
)

# Open the optimization file:
m = read(f"{FILE}")
# name=FILE_NAME.strip('.mps')
# # Define MIP start, if exists:
# FILE_NAME_START = f"{FILE_NAME}.mst"
# if os.path.isfile(FILE_NAME_START):
#     print(">>>> loading {FILE_NAME_START} file")
#     m.read(FILE_NAME_START)

# # Define sol start, if exists:
# FILE_NAME_sol_START = f"{FILE_NAME}.sol"
# if os.path.isfile(FILE_NAME_sol_START):
#     print(">>>> loading {FILE_NAME_sol_START} file")
#     m.read(FILE_NAME_sol_START)

# Define output file names:
log_name = f"{FNAME}.txt"
sol_name = f"{FNAME}.sol"
mps_name = f"{FNAME}_out.mps"

# if write_sols != 0:
#     sol_write = f"{FILE_NAME}_{SEED_OPTIM}"
# else:
# sol_write = ""

# Define settings for the Gurobi optimization:
m.params.LogFile = log_name
m.params.OutputFlag = 1
m.params.TimeLimit = (60 * 60 * HOURS_OPTIM) +60*60  # time limit (s) x min for small jobs
m.params.NonConvex = 2
m.params.Threads = CPUS_OPTIM
# m.params.Seed = SEED_OPTIM
# m.params.SolFiles = sol_write
m.params.FeasibilityTol = 1e-5
m.params.IntFeasTol = 1e-5
m.params.OptimalityTol = 1e-5
m.params.FuncPieces=-2
m.params.FuncPieceError=1e-6
m.params.FuncNonlinear=1
n_scenarios=m.NumScenarios

# def mycallback(model, where):
#     if where == GRB.Callback.MIPSOL:
#         print("solution found", model.cbGet(GRB.Callback.MIPSOL_OBJ))
#         with open(f"res_{FILE_NAME}.txt", "a") as myfile:
#             time = model.cbGet(GRB.Callback.RUNTIME)
#             objval = model.cbGet(GRB.Callback.MIPSOL_OBJ)
#             myfile.write(f"\n {time}_{objval}")


# def callback_VA(model,where):
#     def _callback_save_sol(m):
#         obj_val_dict = {}
#         for i in range(0, m.NumScenarios, 2):
#             rxn_idx = int(i / 2)
#             m.params.ScenarioNumber = i
#             m.update()
#             try :
#                 gur_objval = m.cbGet(GRB.Callback.MIPSOL_OBJ)
#             except:
#                 gur_objval = float('nan')
#             obj_val_dict[rxn_idx] = [(-1) * gur_objval]
#             #max
#             m.params.ScenarioNumber = i + 1
#             m.update()
#             try : 
#                 gur_objval = m.cbGet(GRB.Callback.MIPSOL_OBJ)
#             except:
#                 gur_objval = float('nan')
#             obj_val_dict[rxn_idx].append((+1) * gur_objval)

#         with open(f"{FILE_NAME}_tmp_objval.txt", "w") as f:
#             for k, val in obj_val_dict.items():
#                 f.writelines(f"{k}: {val}\n")


#     if where == GRB.Callback.MIPSOL:
#         # if GRB.Callback.MIPSOL_OPENSCENARIOS < n_scenarios/2:
#             _callback_save_sol(model)




def save_multiscenario_solutions(m):
    #save multiscenario solutions
    no_scenarios = m.NumScenarios
    if no_scenarios > 0:
        obj_val = {}
        obj_bound = {}
        optimal_bounds = {}
        MIPGaps = {}
        for i in range(0, no_scenarios, 2):
            rxn_idx = int(i / 2)
            # Minimization:
            m.params.ScenarioNumber = i
            m.update()
            ObjBound = m.ScenNObjBound
            ObjVal = m.ScenNObjVal
        #  print(rxn_idx, ObjBound, ObjVal)
            if ObjVal != 0:
                MIPGap = abs((ObjBound-ObjVal)/ObjVal)
            else:
                MIPGap = 0

            obj_val[rxn_idx] = [(-1) * ObjVal]
            obj_bound[rxn_idx] = [(-1) * ObjBound]
            MIPGaps[rxn_idx] = [MIPGap]

            if MIPGap <= m.params.MIPGap:
                optimal_bounds[rxn_idx] = [(-1) * ObjBound]
            else:
                optimal_bounds[rxn_idx] = [float('nan')]

            # Maximization:
            m.params.ScenarioNumber = i + 1
            m.update()
            ObjBound = m.ScenNObjBound
            ObjVal = m.ScenNObjVal
            if ObjVal != 0:
                MIPGap = abs((ObjBound-ObjVal)/ObjVal)
            else:
                MIPGap = 0
            obj_val[rxn_idx].append((+1) * m.ScenNObjVal)
            obj_bound[rxn_idx].append((1) * m.ScenNObjBound)
            MIPGaps[rxn_idx].append(MIPGap)
            if MIPGap <= 0.0001:
                optimal_bounds[rxn_idx].append((1) * ObjBound)
            else:
                optimal_bounds[rxn_idx].append(float('nan'))

        with open(f"{OUTPUT_FOLDER}/{FNAME}_objval.txt", "w") as f:
            for k, val in obj_val.items():
                f.writelines(f"{k}: {val}\n")

        # with open(f"{FILE_NAME}_objbounds.txt", "w") as f:
        #     for k, val in obj_bound.items():
        #         f.writelines(f"{k}: {val}\n")

        # with open(f"{FILE_NAME}_optimal_bounds.txt", "w") as f:
        #     for k, val in optimal_bounds.items():
        #         f.writelines(f"{k}: {val}\n")

        # with open(f"{FILE_NAME}_MIPGaps.txt", "w") as f:
        #     for k, val in MIPGaps.items():
        #         f.writelines(f"{k}: {val}\n")

#dont use callbacks for multi scenario as this might slow optimization
if m.NumScenarios > 0:
    #increase MIP gap for multiscenario variability anlaysis as very small suboptimal bounds are probably ok
    m.params.MIPGap = 0.001
    m.optimize()
else:
    m.optimize(mycallback)
m.write(sol_name)

save_multiscenario_solutions(m)
