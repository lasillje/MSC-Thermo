## Folder Structure

- tfva_ecoli.ipynb
    - Main jupyter notebook containing all the code used to generate TVA / FVA / TFVA cluster files, results, and graphs

- setup_VA.sh, gurobi_optimization_VA.py
    - Bash and python script used on the cluster to perform TVA / TFVA

- datafiles/
    - Contains all necessary files to create and load the stoichiometric E. coli model
- results/
    - Contains data from cluster for results of metabolite TVA, and TFVA for Glucose/Acetate growth conditions
    - results/bounds/
        - Contains intermediate results for FVA / TFVA bounds processing for the creation of the composite results
    - results/bounds/composite
        - Contains the flux bounds for the metabolome/no metabolome data and the corresponding composite bounds (*_final_tfva.csv) for both growth conditions
    
