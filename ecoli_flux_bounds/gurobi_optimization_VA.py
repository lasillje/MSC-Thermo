import sys
import os

from gurobipy import read, GRB


print(sys.argv)

# Parse the input arguments:
FILE = str(sys.argv[1])
FNAME=sys.argv[2]
HOURS_OPTIM = int(sys.argv[3])
CPUS_OPTIM = int(sys.argv[4])
OUTPUT_FOLDER=sys.argv[5]
MIPGAP = float(sys.argv[6])

print(
    f"Filename: {FILE}\nDuration of optimization: {HOURS_OPTIM} h\nCPUs used: {CPUS_OPTIM}"
)

m = read(f"{FILE}")

log_name = f"{FNAME}.txt"
sol_name = f"{FNAME}.sol"
mps_name = f"{FNAME}_out.mps"

m.params.LogFile = log_name
m.params.OutputFlag = 1
m.params.TimeLimit = (60 * 60 * HOURS_OPTIM) +60*60  # time limit (s) x min for small jobs
m.params.Threads = CPUS_OPTIM
n_scenarios=m.NumScenarios

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

if m.NumScenarios > 0:
    if m.NumScenarios > 1:
        m.params.MIPGap = MIPGAP
    else:
        m.params.MIPGap = MIPGAP
    m.optimize()
else:
    m.optimize()
    
m.write(sol_name)
save_multiscenario_solutions(m)
