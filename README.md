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

## Step 2 - Run a code to generate mechanistic data
