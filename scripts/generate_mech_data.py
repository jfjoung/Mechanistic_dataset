import os
import logging
import json
import pickle
from tqdm import tqdm
from multiprocessing import Pool
from templates import Reaction_templates
from utils import reaction_process
from utils.elementary_reactions import Get_Reactions
from utils.exceptions import *
from utils.reaction_process import flatten_list
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

def merge_dicts(d, u):
    for k, v in u.items():
        if isinstance(v, dict):
            d[k] = merge_dicts(d.get(k, {}), v)
        else:
            d[k] = d.get(k, 0) + v
    return d

def calling_rxn_template(reaction_dict):
    """
    Retrieve the elementary reaction templates of each reaction
    reaction_dict: one reaction in reactions_with_conditions
    """
    rxn_class_name = reaction_dict['reaction_name']
    rxn_condition = reaction_dict['conditions']
    class_key = get_class_key(rxn_class_name)
    rxn_templates = Reaction_templates.class_reaction_templates[class_key]

    conditions = [rxn_templates[condition] for condition in rxn_condition]

    return conditions

def get_class_key(class_name):
    """
    returns the key(tuple) of class_reaction_templates that includes the class_name
    """

    for classes in Reaction_templates.class_reaction_templates.keys():
        if class_name in classes:
            return classes

    return None

def reagent_matching_for_single_reaction(reaction, class_key):
    '''
    To match reagents for a certain class from list of reaction, in the form of {'rxnsmiles': {'reaction_name': name}}.
    reactions: list of reactions
    '''
    reactants, agents, products = reaction['reaction_smiles'].split(">")
    mols = [Chem.MolFromSmiles(smi, sanitize=False) for smi in reactants.split(".") + agents.split(".")]
    reaction['conditions'] = []
    class_key = get_class_key(class_key)
    if class_key in Reaction_templates.class_reaction_templates:
        class_data = Reaction_templates.class_reaction_templates[class_key]
        for cond_name, cond_data in class_data.items():
            if cond_data['Reagent'] or cond_data['Exclude_reagent']:
                cond_mols = []
                if cond_data['Reagent']:
                    cond_mols = [Chem.MolFromSmarts(x) for x in cond_data['Reagent']]
                exclude_cond_mols = []
                if cond_data['Exclude_reagent']:
                    exclude_cond_mols = [Chem.MolFromSmarts(x) for x in cond_data['Exclude_reagent']]

                matched_reagents = []
                exclude = False

                for patt in cond_mols:
                    patt.UpdatePropertyCache(strict=False)
                    for mol in mols:
                        mol.UpdatePropertyCache(strict=False)
                        #print('Mol :', Chem.MolToSmiles(mol), 'Pat :', Chem.MolToSmarts(patt))
                        if mol.GetSubstructMatch(patt):
                            matched_reagents.append(mol)
                            #print('Matched Mol :', Chem.MolToSmiles(mol), 'Pat :', Chem.MolToSmarts(patt))
                            break
                for patt in exclude_cond_mols:
                    patt.UpdatePropertyCache(strict=False)
                    for mol in mols:
                        mol.UpdatePropertyCache(strict=False)
                        if mol.GetSubstructMatch(patt):
                            #print('Matched Mol :', Chem.MolToSmiles(mol), 'Pat :', Chem.MolToSmarts(patt))
                            #print("Excluding")
                            exclude = True
                            break
                    if exclude: break
                if not exclude and len(cond_mols) == len(matched_reagents):
                    #logging.info(f"Matched: {cond_name}")
                    reaction['conditions'].append(cond_name)
            else:
                reaction['conditions'].append(cond_name)
    return reaction