Fitting of signatures with deconstructSigs v.1.8.0

------------------------Install deconstructSigs with conda ------------------------
Create a conda enviroment to install R and deconstructSigs.
Here are the instructions:
https://github.com/raerose01/deconstructSigs


The script assign_signature_to_mutation.py calls the R script to run deconstructSigs. It also
activates the enviroment where it is installed as <source activate deconstructR>

You can create the environment with deconstructSigs with the deconstructR.yml file

conda env create -f deconstructR.yml

Subsets of signatures to input to the python script are defined as:

BALL_subset: "SBS1","SBS5","SBS18", "SBS2", "SBS13","SBS36","SBS37","SBS9"
TALL_subset: "SBS1","SBS5","SBS18"
ALL_primary_subset: "SBS1", "SBS2","SBS5", "SBS6","SBS9", "SBS13", "SBS17a", "SBS17b","SBS18","SBS34","SBS36", "SBS37"
TALL_HSCP_comparative: "SBS1","SBS5","SBS18", "SBS_hscp"
TALL_relapse_subset: "SBS1","SBS5", "SBS18","SBS32", "SBSA_new", "SBSB_new"

How the inputs were created and examples of how to run it are in make_inputs_fitting_adults.ipynb
