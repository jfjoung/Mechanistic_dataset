import os
import logging
import json
import pickle
from multiprocessing import Pool, Value
import traceback

# Third-party libraries
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdchem, rdmolops
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import networkx as nx
from tqdm import tqdm
from replace import reactions

def compileReactions(rxns): 
    compiledRxns = {}
    for key, value in rxns.items():
        if type(key) is str:
            reactants = (Chem.MolFromSmarts(key),)
        else:
            reactants = tuple(Chem.MolFromSmarts(mol) for mol in key)
        reaction = AllChem.ReactionFromSmarts(value)
        compiledRxns[reactants] = reaction
    
    return compiledRxns

def find_matched_reactions(mols, rxns):
    matched_reactions = {}

    for reactants, reaction in rxns.items():
        used_mols = set()
        found_reactants = []

        for reactant in reactants:
            for mol in mols:
                if mol in used_mols:
                    continue
                if mol.GetSubstructMatch(reactant):
                    used_mols.add(mol)
                    found_reactants.append(mol)
                    # print(f"Mol: {Chem.MolToSmiles(mol)}   Pat: {Chem.MolToSmiles(reactant)}")
        
        # print(found_reactants)
        # print(reactants)
        if len(found_reactants) == len(reactants):
            matched_reactions[reaction] = found_reactants
            #print(matched_reactions)

    return matched_reactions

def preprocess(input, output, rxns):
    """
    Reads a txt file where each line contains:
    <reaction_smiles> <reaction_description>
    """
    with open(input, 'r') as f, open(output, 'a') as o:
        for i, line in enumerate(tqdm(f)):
            line = line.strip()
            if not line or line.startswith('#'):
                continue  # skip empty lines or comments
            # Split on first whitespace or tab
            parts = line.split(maxsplit=1)
            if len(parts) < 2:
                logging.warning(f"Line skipped (not enough parts): {line}")
                continue
            smiles = parts[0]
            r, r1, products = smiles.split(">")
            mols = [Chem.MolFromSmiles(smi, sanitize=False) for smi in r.split(".") + r1.split(".")]
            for mol in mols:
                mol.UpdatePropertyCache(strict=False)
            reactions = find_matched_reactions(mols, rxns)

            if reactions: 
                used_reactants = []
                new_reactants = []
                err = False
                
                for reaction, reactants in reactions.items():
                    try:
                        new_reactants += [Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True) for mol in reaction.RunReactants(reactants)[0]]
                        used_reactants += reactants
                    except Exception as e:
                        err = True  
            
                unused_reactants = [Chem.MolToSmiles(mol) for mol in mols if mol not in used_reactants]
                new_reaction = '.'.join(new_reactants + unused_reactants) + ">>" + products
                o.write(new_reaction.replace("P(", "[P](") + " " + parts[1] + '\n')
                    #o.write(new_reaction + " " + parts[1] + '\n')
                    #print(new_reaction + " " + parts[1] + '\n')

            else:
                o.write(line + "\n")

if __name__ == '__main__':
    #input_file = '../dataset_balancing/h_m_s.txt'
    input_file = '../data/barluenga.txt'
    base, ext = os.path.splitext(input_file)
    output_file = base + '_processed' + ext
    # reactants = preprocess(input_file)
    # print(f"Read {len(reactants)} reactant entries.")
    rxns = compileReactions(reactions)
    preprocess(input_file, output_file, rxns)
