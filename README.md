# A code for Mechanistic Dataset generation
A python code for mechanistic dataset generation. 

## Step 0 - Installation
```shell
git clone https://github.com/jfjoung/Mechanistic_dataset.git
cd Mechanistic_dataset
conda env create --file environment.yml
conda activate mechdata
```

## Step 1 - Data preparation
In order to generate mechanistic dataset from the overall reactions, reaction SMILES and its reaction class name as in NameRXN are required. The reaction SMILES and its reaction class name should be separated by a space as shown below.
It doesn't require atom mapping, a code will automatically assign all the atom mapping number. 
```text
BrCCCOC1=CC(OCO2)=C2C=C1.OC3C(CCN(C3)C(OC(C)(C)C)=O)C4=CC=C(O)C=C4.[OH-]>>OC1C(CCN(C1)C(OC(C)(C)C)=O)C2=CC=C(OCCCOC3=CC(OCO4)=C4C=C3)C=C2 Williamson ether synthesis
```

## Step 3 - Set configuation parameters

Configuration Parameters
This section outlines the configuration parameters used in the script, detailing their purpose and functionality.

Parameters Overview
num_combination: Limits the number of reactant combinations per template to prevent combinatorial explosion.
uni_rxn: Adds a unimolecular version of a bimolecular reaction template while retaining the original.
proton: Ensures proton balance by adjusting templates on-the-fly based on the presence of acids or bases.
max_num_temp: Limits the number of possible templates when modifications like uni_rxn or proton are applied.
stoichiometry: Allows duplicated reactants to be added if stoichiometry errors are detected.
do_not_pruning: Controls whether to prune the reaction network to retain only those that lead to the final product.
num_cycles: Sets a limit on the number of catalytic cycles to avoid long computation times during conversion to SMARTS.
num_reaction_node: Limits the number of reactions in the reaction network to prevent excessive processing time.
byproduct: Adds byproducts from previous reactions when converting to elementary reaction SMARTS.
spectator: Includes all chemical species in the flask at the current step when converting to SMARTS.
full: Determines whether to return overall reaction SMARTS instead of elementary reaction SMARTS.
end: Adds reactions that generate a product from an already existing product.
plain: Generates reaction SMARTS without atom mapping numbers.
explicit_H: Explicitly represents hydrogens in the reaction SMARTS.
reagent: Positions non-reacting species between '>>' in Reaction SMARTS (Reactants > Reagents > Products).
remapping: Remaps atom mapping in SMARTS to start from 1.
data: Path to the file containing overall reactions.
save: Path to the file where results will be saved.
debug: Logs reactions with errors for debugging purposes.
all_info: Saves all reactions, including the reaction network.
stat: Generates statistics for the reactions.
rxn_class: Specifies a particular reaction class to extract elementary reactions from.
process: Number of processors to use in multiprocessing.
verbosity: Sets the logging level.

## Step 3 - Run a code to generate mechanistic data

Mechanistic data can be generated with either of the two methods: \
1. with bash script\
<code>generate_mechanism_data.sh</code> can be modified with the arguments. 
```shell
bash run/generate_mechanism_data.sh
```

2. with python script
```shell
python run.py
```
