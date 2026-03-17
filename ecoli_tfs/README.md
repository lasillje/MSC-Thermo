## Folder Structure

- datafiles/
    - Contains all necessary files to create and load the stoichiometric E. coli model
- results/
    - Contains all intermediate results of TFVA and the TFS files in preperation for the thermodynamically constrained sampling.
    - results/bounds/
        - Final composite TFVA bounds previously obtained from the ecoli_flux_bounds scripts
    - results/tfs
        - Contains necessary files to recreate the thermodynamic-stoichiometric model for the PTA package to run TFS
- escher/
    - Contains files to recreate the Escher map for the acetate growth condition for inspection of Hexokinase (HEX1) activity
- tfs_scripts/
    - Contains modified PTA package python files
    - Introduces the option for linear constraints for solving the PMO problem
    - Allows for manual specification of metabolite concentrations and reaction energies using Thermo-Flux models and data
