import numpy as np
import pandas as pd
import time
import json
import sys
import os
import csv
from multiprocessing import Pool
from tqdm import tqdm

home_dir=os.getcwd()
os.chdir('../..')

import askcos.global_config as gc
from askcos.utilities.canonicalization import canonicalize
from askcos.synthetic.reaction_classification.reaction_class import ReactionClass

os.chdir(home_dir)


def remove_atom_mapping(smiles):
    
    
    if ">" in smiles:
        # Reaction
        try:
            reactants, agents, products = smiles.split(">")
        except ValueError:
            return smiles

        reactants = ".".join(
            sorted(
                remove_atom_mapping(smi)
                for smi in reactants.split(".")
            )
        )
        products = ".".join(
            sorted(
                remove_atom_mapping(smi)
                for smi in products.split(".")
            )
        )
        agents = ".".join(
            sorted(
                remove_atom_mapping(smi)
                for smi in agents.split(".")
            )
        )
        return reactants + ">" + agents + ">" + products
    
    else:
        # Molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(0)
            return Chem.MolToSmiles(mol)
        
        return smiles

    
    

def reaction_classification(rxn_list, reaction_classifier):
    
    
    results = reaction_classifier.predict(rxn_list, num_results=1)
    rxn_dicts = []
    for rxnsmi, result in zip(rxn_list, results): 
        rxn_dicts.append({rxnsmi: result['result'][0]})
    
    # rxn_dicts: list of {rxnsmiles: rxn_class_data}
    # rxn_class_data ex: {'rank': 1, 
    #                     'reaction_num': '3.1.3', 
    #                     'reaction_name': 'Iodo Suzuki coupling', 
    #                     'reaction_classnum': '3.1', 
    #                     'reaction_classname': 'Suzuki coupling', 
    #                     'reaction_superclassnum': '3', 'reaction_superclassname': 
    #                     'C-C bond formation', 
    #                     'prediction_certainty': 0.8232529759407043}}

    return rxn_dicts

def batched_reaction_classification(rxn_list, batch_size = 1000):
    
    #p = Pool(5)
    batches = [rxn_list[i:i+batch_size] for i in range(0, len(rxn_list), batch_size)]
    
    reaction_classifier = ReactionClass()
    reaction_classifier.load_model()
    
    classified_rxns = []
#     for result in tqdm(p.imap(reaction_classification, batches), total=len(batches)):
#         classified_rxns = classified_rxns+result
    
#     p.close()
#     p.join()

    logfile = f'../data/classified_USPTO_{datatype}_log.json'
    
    for i, rxn_batch in enumerate(batches):
        classified_rxns = classified_rxns + reaction_classification(rxn_batch, reaction_classifier)
    
        print(f"{i}/{len(batches)}")
        with open(logfile, 'w') as out_file:
            json.dump(classified_rxns, out_file)
    
    
    return classified_rxns
    

if __name__ == '__main__':
    
    datatype = 'MIT'
    file_path = f'../data/USPTO_{datatype}.csv'
    data = pd.read_csv(file_path)
    
    start = time.time()
    
    print("Classification Started")

    result = batched_reaction_classification(list(data['reactions']))
    print("Elapsed Time:", time.time() - start)
    print("Num. Reactions", len(result))
    output_file = f'../data/classified_USPTO_{datatype}.json'
    with open(output_file, 'w') as out_file:
        json.dump(result, out_file)