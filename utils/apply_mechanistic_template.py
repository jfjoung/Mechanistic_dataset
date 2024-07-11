import os
import sys
import json
import rdkit
import gzip
import logging
import copy
import numpy as np
from . import AcidBase_lookup, Reaction_templates
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
from joblib import Parallel, delayed
import itertools
from itertools import product, permutations
import networkx as nx
from collections import defaultdict
import matplotlib.pyplot as plt
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def mols_from_smiles_list(all_smiles):
    '''Given a list of smiles strings, this function creates rdkit
    molecules'''
        
    mols = []
    for smiles in all_smiles:
        if not smiles: continue
        mols.append(Chem.MolFromSmiles(smiles,sanitize=False))
    return mols

def remove_atom_map(mol, isotope=False):
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
        if isotope: atom.SetIsotope(0)
    return mol

def isotope_to_atommap(mol):
    for idx, a in enumerate(mol.GetAtoms()):
        a.SetAtomMapNum(a.GetIsotope())
        a.SetIsotope(0)
    return mol

def get_class_key(class_name):

    '''
    returns the key(tuple) of class_reaction_templates that includes the class_name
    '''

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
    mols = [Chem.MolFromSmiles(smi) for smi in reactants.split(".") + agents.split(".")]
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
                    for mol in mols:
                        if mol.GetSubstructMatch(patt):
                            matched_reagents.append(mol)
                            break
                for patt in exclude_cond_mols:
                    for mol in mols:
                        if mol.GetSubstructMatch(patt):
                            exclude = True
                            break
                    if exclude: break
                if not exclude and len(cond_mols) == len(matched_reagents):
                    reaction['conditions'].append(cond_name)
            else:
                reaction['conditions'].append(cond_name)
    return reaction

def calling_rxn_template(reaction_dict):
    """
    Retrive the elementary reaction templates of each reaction
    reaction_dict: one reaction in reactions_with_conditions
    """
    rxn_class_name=reaction_dict['reaction_name']
    rxn_condition=reaction_dict['conditions']
    class_key=get_class_key(rxn_class_name)
    rxn_templates=Reaction_templates.class_reaction_templates[class_key]

    conditions=[rxn_templates[condition] for condition in rxn_condition]

    return conditions

def prepare_reactants(reaction_dict, args):
    """
    Preprocessing to make reactants with isotope labeling
    reaction_dict: one reaction in reactions_with_conditions
    """

    rxn_smi=reaction_dict['reaction_smiles']

    reactants, agents, products = [mols_from_smiles_list(x) for x in
                                [mols.split('.') for mols in rxn_smi.split('>')]]

    rmol=reactants+agents

    idx=1
    reactant_dict=dict()
    _smiles_list=[]

    # if args.explicit_H:
    #     ps = Chem.SmilesParserParams()
    #     ps.removeHs = False
    # else:
    ps = Chem.SmilesParserParams()


    for num, mol in enumerate(rmol):
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol), ps)
        mol = remove_atom_map(mol)
        # if args.explicit_H:
        # mol = Chem.AddHs(mol, explicitOnly=False)
        start_idx=idx
        for atom in mol.GetAtoms():
            atom.SetIsotope(idx)
            idx+=1

        _smiles_list.append(Chem.MolToSmiles(mol))

        #Save molecule in a dictionary

        reactant_dict[num] = {'smiles_w_isotope':Chem.MolToSmiles(mol),  # Using isotope as atom mapping for reaction
                            'atom_mapping':[i for i in range(start_idx,idx)],  #Checking for atom-mapping collision
                           'smiles_w_mapping': Chem.MolToSmiles(isotope_to_atommap(mol)),  # Output for generating elementary reaction graph
                              'smiles': Chem.MolToSmiles(remove_atom_map(mol, isotope=True)),  # Plain SMILES string
                           'rxn_history': list(),
                            'identity': 'reactant'
                              }

    reactant_dict['last_mapping_number']=idx-1
    return reactant_dict

def check_products_validity(outcome):
 #Many outcomes can be produced from one template (e.g. symmetric molecules)

    good_molecule=0 # To check every molecule is making sense
    for prod_mol in outcome:
        try:# Not realistic molecule should be rejected, except aromaticity error.
            prod_mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(prod_mol,
                             Chem.SanitizeFlags.SANITIZE_FINDRADICALS|Chem.SanitizeFlags.SANITIZE_SETAROMATICITY|Chem.SanitizeFlags.SANITIZE_SETCONJUGATION|Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION|Chem.SanitizeFlags.SANITIZE_SYMMRINGS,catchErrors=True)
            good_molecule+=1

        except:
            continue
    if len(outcome)==good_molecule:
        return True
    else: return False

# def remove_isotope(prod_mol):
#     prod_istope_smi=Chem.MolToSmiles(prod_mol)
#     for a in prod_mol.GetAtoms():
#         a.SetIsotope(0)
#     return Chem.MolFromSmiles(prod_istope_smi), Chem.MolToSmiles(prod_mol)

def product_dict(prod_mol, reactant_id, reactant_history, args):
    prod_smi_isotope = Chem.MolToSmiles(prod_mol) # Get isotompe-mapped smiles
    atom_mapping=list()

    # if args.explicit_H:
    #     ps = Chem.SmilesParserParams()
    #     ps.removeHs = False
    # else:
    ps = Chem.SmilesParserParams()

    for a in prod_mol.GetAtoms(): #Change isotope to atom map
        a.SetAtomMapNum(a.GetIsotope())
        atom_mapping.append(a.GetIsotope())
        a.SetIsotope(0)
    prod_smi_map=Chem.MolToSmiles(prod_mol) #Get atom mapped smiles
    for a in prod_mol.GetAtoms(): #Remove any mapping
        a.SetAtomMapNum(0)

    # Get plain smiles
    try:
        prod_smi = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(prod_mol), ps), isomericSmiles=False)
    except:
        prod_smi = Chem.MolToSmiles(prod_mol, isomericSmiles=False)

    reactant_history = list(set([item for sublist in reactant_history if isinstance(sublist, set) for item in sublist]+\
                        [item for item in reactant_history if isinstance(item, int)]+\
                        [item for sublist in reactant_history if isinstance(sublist, list) for item in sublist]))

    new_prod={'smiles_w_isotope':prod_smi_isotope,
            'atom_mapping':atom_mapping,
            'smiles':prod_smi,
           'smiles_w_mapping': prod_smi_map,
           'rxn_history': reactant_history,
            'identity': 'intermediate'
             }
    return new_prod

def proton_balanced_template(rxn_flask, pKas, rxn_templates, args):

    """
    rxn_flask: dictionary - key: species index (0~), value: dictionary of species data (smiles, tags, etc)
    pKas: list of pKas for each rxn_template in rxn_templates
    """
    proton_balanced_rxn_templates = []
    original_rxn_templates = []

    if pKas == [None]*len(pKas):
#         for data in rxn_flask.values():
#             denotes whether a species works as acid or base in each template
#             data['acid_or_base'] = [None for temp in rxn_templates]
#         return rxn_flask, rxn_templates
        return rxn_templates, proton_balanced_rxn_templates

    for pKa, templ in zip(pKas, rxn_templates):
        if not pKa:
            original_rxn_templates.append(templ)
            continue
        A=pKa.get('A')
        B=pKa.get('B')
        possible_acid_base=[]
        if A:
            filtered_data = [d for d in AcidBase_lookup.Acid_base if 'A' in d['role'] and d['pKa'] <= A]
            sorted_data = sorted(filtered_data, key=lambda x: x['pKa'])
            possible_acid_base=find_acid_base(rxn_flask, sorted_data, 'A', args)
        if B:
            filtered_data = [d for d in AcidBase_lookup.Acid_base if 'B' in d['role'] and d['pKa'] >= B]
            sorted_data = sorted(filtered_data, key=lambda x: x['pKa'], reverse=True)
            possible_acid_base=find_acid_base(rxn_flask, sorted_data, 'B', args)
        if not possible_acid_base: continue

        for acid_base in possible_acid_base:
            # for templ in rxn_templates:
            r, p = templ.split('>>')
            new_r = '.'.join([r, acid_base[0]])
            new_p = '.'.join([p, acid_base[1]])
            new_templ='>>'.join([new_r,new_p])
            proton_balanced_rxn_templates.append(new_templ)

    return original_rxn_templates, proton_balanced_rxn_templates

def allow_unimolecular_rxn(rxn_templates):
    new_templates=[]
    for templ in rxn_templates:
        rxn = AllChem.ReactionFromSmarts(templ)
        num_reactants=rxn.GetNumReactantTemplates()

        # new_templates.append(templ)
        if num_reactants > 1:
            r, p = templ.split('>>')
            r= '('+r+')'
            p='('+p+')'
            new_temp='>>'.join([r,p])
            new_templates.append(new_temp)
    return new_templates

def find_acid_base(rxn_flask, filtered_list, ab_condition, args):
    possible_acid_base=[]

    # if args.explicit_H:
    #     ps = Chem.SmilesParserParams()
    #     ps.removeHs = False
    #     ps.sanitize = False
    # else:
    ps = Chem.SmilesParserParams()
    ps.sanitize = False

    for acid_base in filtered_list:
        if ab_condition=='A':
            reactant=acid_base['Acid']
            product=acid_base['Base']
        elif ab_condition=='B':
            reactant=acid_base['Base']
            product=acid_base['Acid']

        patt = Chem.MolFromSmarts(reactant)
        mols = [Chem.MolFromSmiles(rxn_flask[x]['smiles'], ps) for x in rxn_flask  if type(x) is int]
        for mol in mols:
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(mol,
                             Chem.SanitizeFlags.SANITIZE_FINDRADICALS | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                             catchErrors=True)
            if mol and mol.GetSubstructMatch(patt) and [reactant, product] not in possible_acid_base:
                possible_acid_base.append([reactant, product])
    return possible_acid_base

def find_reactants(rxn_flask, rxn_templates, args):

    # if args.explicit_H:
    #     ps = Chem.SmilesParserParams()
    #     ps.removeHs = False
    #     ps.sanitize = False
    # else:
    ps = Chem.SmilesParserParams()
    ps.sanitize = False

    reactive_dict = dict()
    mol_dict = {reactant: Chem.MolFromSmarts(rxn_flask[reactant]['smiles_w_mapping']) for reactant in rxn_flask if type(reactant) is int}

    for num_temp, templ in enumerate(rxn_templates):
        reactant_dict=dict()
        rxn = AllChem.ReactionFromSmarts(templ)
        r = [rmol for rmol in rxn.GetReactants()]
        num_reactants=rxn.GetNumReactantTemplates()

        templ_dict = {templ: r}
        # Find reactive chemical species
        for mol_key in mol_dict:
            mol = mol_dict[mol_key]
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(mol,
                             Chem.SanitizeFlags.SANITIZE_FINDRADICALS | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                             catchErrors=True)

            for temp_key in templ_dict:
                patterns = templ_dict[temp_key]
                for pat in patterns:
                    if mol and mol.GetSubstructMatch(pat):
                        if Chem.MolToSmarts(pat) not in reactant_dict:
                            reactant_dict[Chem.MolToSmarts(pat)] = []
                        reactant_dict[Chem.MolToSmarts(pat)].append(mol_key)
        combinations = [list(combination) for combination in itertools.product(*reactant_dict.values()) if len(combination)==num_reactants]
        # Check the atom map collision
        reactive_combination = list()
        reactant_history_per_com=list()
        if args.stoichiometry:
            # To save every combination
            duplicated_combi=list()
        for combination in combinations:
            if num_reactants>1:
                reactant_atommap=[set(i) for i in [rxn_flask[mol_num]['atom_mapping'] for mol_num in combination]]
                # Check if there is the same atom in two reactants
                if reactant_atommap[0]&reactant_atommap[1]:
                    duplicated_combi.append(combination)
                    continue
            reactant_history=[i for i in [rxn_flask[mol_num]['rxn_history'] for mol_num in combination]]
            reactant_history.extend(combination)
            # Check if at least one species came from the same reactant
            if has_duplicates(reactant_history):
                if args.stoichiometry:
                    duplicated_combi.append(combination)
                    continue
                else:
                    continue

            reactive_combination.append(combination)
            reactant_history_per_com.append(reactant_history)

            if num_reactants>1:
                reactive_combination.append(list(reversed(combination)))
                reactant_history_per_com.append(reactant_history)
            reactive_dict[templ]={'combination': reactive_combination,
                              'num_reactants': num_reactants,
                              'reactant_history': reactant_history_per_com}
        '''
        Hereafter, it is a code for matching the stoichiometric details
        '''
        if args.stoichiometry and duplicated_combi and not reactive_dict:
            idx = rxn_flask['last_mapping_number']+1
            new_molecule_dict = {}
            new_combination_list = []
            for combination in duplicated_combi:
                new_combination = []
                for mol_num in combination:
                    if rxn_flask[mol_num]['identity']=='reactant' and mol_num not in new_molecule_dict.keys():
                        start_idx=idx
                        new_mol = Chem.MolFromSmiles(rxn_flask[mol_num]['smiles'], ps)
                        for atom in new_mol.GetAtoms():
                            atom.SetIsotope(idx)
                            idx+=1

                        rxn_flask[len(rxn_flask)+1]={'smiles_w_isotope':Chem.MolToSmiles(new_mol),   # Using isotope as atom mapping for reaction
                            'atom_mapping':[i for i in range(start_idx,idx)],  #Checking for atom-mapping collision
                            'smiles':Chem.MolToSmiles(remove_atom_map(new_mol), isomericSmiles=False),  #Plain SMILES string
                           'smiles_w_mapping': Chem.MolToSmiles(isotope_to_atommap(new_mol), isomericSmiles=False), # Output for generating elementary reaction graph
                           'rxn_history': list(),
                            'identity': 'reactant'
                             }
                        new_molecule_dict[mol_num] = len(rxn_flask)
                        new_combination.append(len(rxn_flask))
                    elif mol_num in new_molecule_dict.keys():
                        new_combination.append(new_molecule_dict[mol_num])
                    else: new_combination.append(mol_num)
                new_combination_list.append(new_combination)

            rxn_flask['last_mapping_number']=idx-1
            for combination in new_combination_list:
                if num_reactants>1:
                    reactant_atommap=[set(i) for i in [rxn_flask[mol_num]['atom_mapping'] for mol_num in combination]]
                    # Check if there is the same atom in two reactants
                    if reactant_atommap[0]&reactant_atommap[1]:
                        continue
                reactant_history=[i for i in [rxn_flask[mol_num]['rxn_history'] for mol_num in combination]]
                reactant_history.extend(combination)
                reactant_history = flatten_list(reactant_history)

                # Check if at least one species came from the same reactant
                if has_duplicates(reactant_history):
                    continue
                reactive_combination.append(combination)
                reactant_history_per_com.append(reactant_history)

                if num_reactants>1:
                    reactive_combination.append(list(reversed(combination)))
                    reactant_history_per_com.append(reactant_history)
            reactive_dict[templ] = {'combination': reactive_combination,
                                    'num_reactants': num_reactants,
                                    'reactant_history': reactant_history_per_com}

    final_reactive_dict=dict()
    for key, value in reactive_dict.items():
        num_reactive_comb=len(value['combination'])
        if num_reactive_comb > args.num_combination:
            if args.verbosity > 2:
                logging.info('The template below has {} reactive combinations  in this reaction flask'.format(num_reactive_comb))
                logging.info('{}'.format(key))
                logging.info('This template is abandoned to prevent combination explosion.')
        else:
            final_reactive_dict[key]=value
    if args.verbosity > 1 and not final_reactive_dict:
        logging.info('No reactants are found!')

    return rxn_flask, final_reactive_dict

def has_duplicates(input_list):
    lst=[]
    elem=[]
    for item in input_list:
        if type(item) is list:
            lst.append(item)
        else:
            elem.append(item)
    for item in elem:
        for lst_item in lst:
            if item in lst_item:
                return True
    return False

def run_single_reaction(rxn_flask, single_step, args):
    """Run single elementary reactions
    rxn_flask: a dictionary contains all the reactants
    single_step: a dictionary contains the templates, pKa condition, and description (elementary step)
    proton: boolean to balance proton
    uni_rxn: boolean to allow unimolecular reaction
    stoichiometry: boolean to allow duplicate reactants
    v: verbose
    """
    rxn_templates = single_step['Templates']
    pKas = single_step['pKa']
    # print("TEMPLATES", rxn_templates)
    reaction_network=[]
    reaction_pair=[]
    if args.verbosity > 2:
        num_original_temp=len(rxn_templates)
        logging.info('The number of the original templates is {}'.format(num_original_temp))

    if args.proton: # Get proton balanced reaction template
        # TODO: label acids and bases (e.g. 'acid' instead of just 'reactant' when labelling), add isAcid, isBase
        rxn_templates, rxn_proton_templates = proton_balanced_template(rxn_flask, pKas, rxn_templates, args)
        num_proton_temp = len(rxn_proton_templates)
        if args.verbosity > 2:
            logging.info('New proton balanced {} templates are generated'.format(max(num_proton_temp, 0)))

    if args.uni_rxn: # Get unimolecular reaction template, Intramolecular proton transfer is not allowed
        rxn_uni_templates = allow_unimolecular_rxn(rxn_templates)
        num_uni_temp = len(rxn_uni_templates)
        if args.verbosity > 2:
            logging.info('New unimolecular {} templates are generated'.format(max(num_uni_temp,0)))

    # if args.explicit_H:
    #     ps = Chem.SmilesParserParams()
    #     ps.removeHs = False
    #     ps.sanitize = False
    # else:
    ps = Chem.SmilesParserParams()
    ps.sanitize = False
    if args.proton:
        rxn_templates = rxn_templates + rxn_proton_templates
    if args.uni_rxn:
        rxn_templates = rxn_templates + rxn_uni_templates

    if len(rxn_templates) > args.max_num_temp:
        if args.verbosity > 2:
                logging.info('{} templates are generated, skipping it to avoid combinatorial explosion'.format(len(rxn_templates)))
        raise ValueError("Too many templates are generated.")

    #Define the reactive chemicals
    #rxn_flask may be updated if stoichiometry is on
    rxn_flask, reactive_dict=find_reactants(rxn_flask, rxn_templates, args)
    # print("RXN_FLASK", rxn_flask)
    # print(reactive_dict)

    prod_smi_list = []
    for templ in reactive_dict:
        # print("TEMPL", templ)
        if args.verbosity > 3:
            logging.info('Template is {}'.format(templ))

        combinations=reactive_dict[templ]['combination']
        num_reactants=reactive_dict[templ]['num_reactants']
        reactant_history=reactive_dict[templ]['reactant_history']
        rxn = AllChem.ReactionFromSmarts(templ)

        for combination, history in zip(combinations, reactant_history):
            reactants=[Chem.MolFromSmiles(rxn_flask[num]['smiles_w_isotope'],ps)
                                 for num in combination]

            outcomes = rxn.RunReactants(reactants)
            # print([rxn_flask[num]['smiles_w_isotope'] for num in combination])
            # print([Chem.MolToSmiles(mol) for mol in reactants])

            # for a in reactants[1].GetAtoms():
            #     print(a.GetIsotope(), a.GetIsAromatic(), a.GetAtomMapNum())
            if not outcomes:
                continue

            if args.verbosity > 3:
                reactants_smi=[rxn_flask[num]['smiles'] for num in combination]
                logging.info('The reactants are {}'.format(reactants_smi))

            for outcome in outcomes:
                reactant_id=set(num for num in combination)
                if not check_products_validity(outcome): continue
                product_id=set()
                if args.verbosity > 3: logging.info('{} products are formed'.format(len(outcome)))

                prod_smi = []
                new_prod_num = []
                for i, prod_mol in enumerate(outcome):
                    new_prod=product_dict(prod_mol, reactant_id, history, args)
                    r_map = set(flatten_list([rxn_flask[rid]['atom_mapping'] for rid in reactant_id]))
                    # prod_nums = [key for key, value in rxn_flask.items() if type(key) is int and value['smiles'] == new_prod['smiles'] and sorted(value['atom_mapping']) == sorted(new_prod['atom_mapping'])]
                    prod_nums = [key for key, value in rxn_flask.items() if
                                type(key) is int and value['smiles_w_mapping'] == new_prod['smiles_w_mapping']]
                    found = False
                    for prod_num in prod_nums:
                        p_map = set(flatten_list([rxn_flask[prod_num]['atom_mapping']]))
                        if prod_num and r_map.issuperset(p_map):
                            found = True
                            if args.verbosity > 3: logging.info('{} is found'.format(rxn_flask[prod_num]['smiles_w_mapping']))
                            product_id.add(prod_num)
                            prod_smi.append(rxn_flask[prod_num]['smiles'])
                            break
                    if not found:
                        num_mols=len(rxn_flask)+1
                        product_id.add(num_mols)
                        new_prod_num.append(num_mols)
                        rxn_flask[num_mols]= new_prod
                        prod_smi.append(rxn_flask[num_mols]['smiles'])
                        if args.verbosity > 3: logging.info('{} is formed'.format(rxn_flask[num_mols]['smiles_w_mapping']))

                prod_smi = set(prod_smi)
                if [reactant_id,product_id] not in reaction_pair and prod_smi not in prod_smi_list:
                    reaction_pair.append([reactant_id,product_id])
                    reaction_network.append([reactant_id,templ, product_id])
                    prod_smi_list.append(prod_smi)
                else:
                    for num in product_id:
                        if num in new_prod_num:
                            del rxn_flask[num]

    return rxn_flask, reaction_network

def run_full_reaction(rxn_flask, condition, cond_name, args):

    """
    Given an initial rxn_flask, applies the templates for each condition and
    returns lists of final reaction flasks and networks for each condition
    Returns empty lists if there are no conditions

    rxn_flasks: a list of final reaction flasks for each condition
    tot_networks: a list of tot_net for each condition
    """

    tot_networks = list()
    rxn_flasks = list()
    for cond, cname in zip(condition, cond_name):
        if args.verbosity > 0:
            logging.info('Start reaction under the condition of {}'.format(cname))
        tot_net = list()
        rxn_flask_cond = copy.deepcopy(rxn_flask)

        for i, temp in enumerate(cond['Stages'].values()):
            if args.verbosity > 1:
                logging.info('Apply {}th template'.format(i))
            rxn_flask_cond, reaction_network=run_single_reaction(rxn_flask_cond, temp, args)
            tot_net+=reaction_network
        tot_networks.append(tot_net)
        rxn_flasks.append(rxn_flask_cond)

    return rxn_flasks, tot_networks


def find_product(example_rxn, rxn_flask, args):

    # if args.explicit_H:
    #     ps = Chem.SmilesParserParams()
    #     ps.removeHs = False
    # else:
    ps = Chem.SmilesParserParams()


    reaction_smi = example_rxn['reaction_smiles']
    _, _, product_smi = reaction_smi.split('>')
    product_smi = Chem.MolToSmiles(Chem.MolFromSmiles(product_smi, ps), isomericSmiles=False)

    product_smi_list=product_smi.split('.')
    real_product_smi_list=[]
    for psmi in product_smi_list:
        pmol=Chem.MolFromSmiles(psmi, ps)
        pat = Chem.MolFromSmarts("[#6]")
        if len(pmol.GetSubstructMatches(pat)) > 0:
            pmol=remove_atom_map(pmol)
            real_product_smi_list.append(Chem.MolToSmiles(pmol))
    # print("Real product SMILES list: ", real_product_smi_list)
    # print("Reaction flask: ", rxn_flask)
    # for key, value in rxn_flask.items():
        # if type(value)==dict: print(value['smiles'])
    matching_keys = [key for key, value in rxn_flask.items() if type(key) is int and value['smiles'] in real_product_smi_list]
    # print("Matching keys: ", matching_keys)
    if matching_keys:
        for key in matching_keys:
            rxn_flask[key]['identity'] = 'product'
        return rxn_flask
    else:
        return None

def topo_pos(G):
    """Display in topological order, with simple offsetting for legibility"""
    pos_dict = {}
    for i, node_list in enumerate(nx.topological_generations(G)):
        x_offset = len(node_list) / 2
        y_offset = 0.0
        for j, name in enumerate(node_list):
            pos_dict[name] = (j - x_offset, -i + j * y_offset)

    return pos_dict

def reaction_network(rxn_flask, tot_network, args):
    '''
    simple: find all shortest paths or all simple paths
    '''
    #Get all the reactants from rxn_flask
    reactant_list=list()
    for key, value in rxn_flask.items():
        if type(key) is int and value['identity']=='reactant':
            reactant_list.append(key)

    # Create a directed graph
    G = nx.DiGraph()
    # Make reaction network
    for num, reaction in enumerate(tot_network):
        first_set = reaction[0]
        template = reaction[1]
        last_set = reaction[-1]

        reaction_node = f"Reaction {num}"
        G.add_node(reaction_node, reaction={'Template': template})
        for chemical in first_set:
            chemical_node=f"Molecule {chemical}"
            molecule=rxn_flask[chemical]
            G.add_node(chemical_node, molecule=molecule)
            G.add_edge(chemical_node, reaction_node)
        for chemical in last_set:
            chemical_node=f"Molecule {chemical}"
            molecule=rxn_flask[chemical]
            G.add_node(chemical_node, molecule=molecule)
            G.add_edge(reaction_node,chemical_node)

    # Process the network
    reactant_ids=list()
    product_ids=list()
    for nid, attrs in G.nodes.data():
        if nid.startswith('Molecule'):
            if attrs['molecule']['identity']=='reactant':
                reactant_ids.append(nid)
            elif not [n for n in G.neighbors(nid)] and attrs['molecule']['identity']!='product':
                attrs['molecule']['identity'] = 'byproduct'
            elif attrs['molecule']['identity']=='product':
                product_ids.append(nid)
    # print(nx.node_link_data(G))
    # Find all the routes connecting the reactants and products
    reaction_path=list()
    for rid in reactant_ids:
        for pid in product_ids:
            try:
                if args.simple:
                    if nx.all_simple_paths(G, source=rid, target=pid):
                        for path in nx.all_simple_paths(G, source=rid, target=pid):
                            if len(reaction_nodes_in_path) < 13:
                                reaction_path.append(node for node in path if node.startswith('Reaction'))
                else:
                    if nx.all_shortest_paths(G, source=rid, target=pid):
                        for path in nx.all_shortest_paths(G, source=rid, target=pid):
                            reaction_nodes_in_path=[node for node in path if node.startswith('Reaction')]
                            if len(reaction_nodes_in_path) < 13:
                                reaction_path.append([node for node in path if node.startswith('Reaction')])
            except: continue

    
    reaction_path = [element for sublist in reaction_path for element in sublist]
    reaction_path = list(set(reaction_path))
    # Get all neighbors of the reaction nodes
    path_node=[]
    for reaction_node in reaction_path:
        path_node.append(reaction_node)
        neighbor_node=[n for n in nx.all_neighbors(G, reaction_node)]
        for nn in neighbor_node:
            path_node.append(nn)

    # TODO: Check if byproducts react further
    if args.do_not_pruning:
        return G
    # Make subgraph of a true path
    G_sub = nx.DiGraph(G.subgraph(path_node))
    # print(nx.node_link_data(G_sub))

    #Check if there is not-connected intermediates.
    not_connected_inter_nodes=list()
    for nid, attrs in G_sub.nodes.data():
        if nid.startswith('Molecule') and attrs['molecule']['identity']== 'intermediate':
            parents=[node for node in G_sub.predecessors(nid)]
            child=[node for node in G_sub.successors(nid)]
            if not parents or not child:
                not_connected_inter_nodes.append(nid)
    # print(nx.node_link_data(G_sub))
    # print('not_connected_inter_nodes', not_connected_inter_nodes)
    if len(not_connected_inter_nodes) > 1:
        for rid, pid in permutations(not_connected_inter_nodes, 2):
            try:
                shortest_paths = list(nx.all_shortest_paths(G, source=rid, target=pid))
                if shortest_paths:
                    new_path = []
                    for path in shortest_paths:
                        checking_path = [node for node in path if node.startswith('Reaction')]
                        if not all(elem in path_node for elem in checking_path):
                            new_path.append(checking_path)
                    new_path=flatten_list(new_path)
                    for reaction_node in new_path:
                        path_node.append(reaction_node)
                        neighbor_node=[n for n in nx.all_neighbors(G, reaction_node)]
                        for nn in neighbor_node:
                            path_node.append(nn)
                    path_node=list(set(path_node))
            except nx.NetworkXNoPath: pass
            try:
                child_rxn_node = set([node for node in G.successors(rid)])
                child_rxn_node2 = set([node for node in G.successors(pid)])
                shared_rxn_node = list(child_rxn_node&child_rxn_node2)
                for rxn_id in shared_rxn_node:
                    child = [node for node in G.successors(rxn_id)]
                    byprode_child = [node for node in G.successors(rxn_id) if G.nodes[node]['molecule']['identity'] == 'byproduct']

                    if set(child) == set(byprode_child) and rxn_id not in path_node:
                        child.append(rid)
                        child.append(pid)
                        child.append(rxn_id)
                        for nn in child:
                            path_node.append(nn)
                path_node = list(set(path_node))
            except: continue
    if not_connected_inter_nodes:
        for int_id in not_connected_inter_nodes:
            rxn_node = set([node for node in G.successors(int_id)])
            for rxn_id in rxn_node:
                child = [node for node in G.successors(rxn_id)]
                byprode_child = [node for node in G.successors(rxn_id) if G.nodes[node]['molecule']['identity'] == 'byproduct']
                non_byprod_child = [node for node in G.successors(rxn_id) if G.nodes[node]['molecule']['identity'] != 'byproduct' and node in list(G_sub.nodes)]
                if set(child) == set(byprode_child+non_byprod_child) and rxn_id not in path_node:
                    child.append(rxn_id)
                    child.append(int_id)
                    for nn in child:
                        path_node.append(nn)

    G_sub = nx.DiGraph(G.subgraph(path_node))

    # print(nx.node_link_data(G_sub))

    not_connected_inter_nodes = True
    i = 0
    while not_connected_inter_nodes:
        not_connected_inter_nodes = []
        for nid, attrs in G_sub.nodes.data():
            if nid.startswith('Molecule') and attrs['molecule']['identity'] == 'intermediate':
                child = [node for node in G_sub.successors(nid)]
                if not child:
                    not_connected_inter_nodes.append(nid)
        # print(not_connected_inter_nodes)
        for nid in not_connected_inter_nodes:
            G_sub.remove_nodes_from([n for n in G_sub.predecessors(nid)])
            G_sub.remove_node(nid)
        i += 1
        if i > args.num_reaction_node:
            break #TODO: Make it raise error.

    not_connected_byprod = True
    while not_connected_byprod:
        not_connected_byprod = []
        for nid, attrs in G_sub.nodes.data():
            if nid.startswith('Molecule') and attrs['molecule']['identity'] == 'byproduct':
                parent = [node for node in G_sub.predecessors(nid)]
                if not parent:
                    not_connected_byprod.append(nid)
            if nid.startswith('Molecule') and attrs['molecule']['identity'] == 'reactant':
                parent = [node for node in G_sub.successors(nid)]
                if not parent:
                    not_connected_byprod.append(nid)
        for nid in not_connected_byprod:
            G_sub.remove_node(nid)
        i += 1
        if i > args.num_reaction_node:
            break #TODO: Make it raise error.

    # print(nx.node_link_data(G_sub))
    # Add spectators
    missing_molecule_nodes = [node_id for node_id in reactant_list if f'Molecule {node_id}' not in G_sub.nodes]

    for node_id in missing_molecule_nodes:
        chemical_node=f"Molecule {node_id}"
        molecule=rxn_flask[node_id]
        molecule['identity']='spectator'
        G_sub.add_node(chemical_node, molecule=molecule)
    # print(nx.node_link_data(G_sub))
    return G_sub

def draw_reaction_graph(G, size=[5,5], labels = True):
    try:
        posit =  topo_pos(G)
    except:
        posit =  nx.kamada_kawai_layout(G, scale=2)
        #pos = nx.nx_agraph.graphviz_layout(G, prog=“dot”)

    plt.figure(3,figsize=(size[0],size[1]))

    mol_list=[]
    reaction_list=[]
    reactant_list=[]
    product_list=[]
    byproduct_list=[]
    spectator_list=[]
    for nid, attrs in G.nodes.data():
        if nid.startswith('Molecule'):
            if attrs['molecule']['identity']=='reactant':
                reactant_list.append(nid)
            elif attrs['molecule']['identity']=='product':
                product_list.append(nid)
            elif attrs['molecule']['identity']=='byproduct':
                byproduct_list.append(nid)
            elif attrs['molecule']['identity']=='spectator':
                spectator_list.append(nid)
            else: mol_list.append(nid)
        else: reaction_list.append(nid)


    nx.draw(G , posit, with_labels = labels)
    nx.draw_networkx_nodes(G, posit, nodelist=mol_list, node_color="tab:red")
    nx.draw_networkx_nodes(G, posit, nodelist=reaction_list, node_shape="s",node_color="tab:green")
    nx.draw_networkx_nodes(G, posit, nodelist=reactant_list,  node_color="tab:blue")
    nx.draw_networkx_nodes(G, posit, nodelist=byproduct_list, node_color="tab:cyan")
    nx.draw_networkx_nodes(G, posit, nodelist=product_list, node_color="magenta")
    nx.draw_networkx_nodes(G, posit, nodelist=spectator_list, node_color="tab:purple")

def get_mechanistic_network(rxn, args):
    condition = calling_rxn_template(rxn)
    if args.verbosity > 2:
        logging.info('{} conditions for this reaction are retrieved'.format(len(condition)))
    rxn_flask=prepare_reactants(rxn, args)
    if args.verbosity > 0:
        logging.info('Start applying the templates')
    rxn_flasks, tot_networks=run_full_reaction(rxn_flask, condition, rxn['conditions'], args)

    G_dict=dict()
    for rxn_condition, rxn_flask,tot_network in zip(rxn['conditions'], rxn_flasks, tot_networks):
        rxn_flask=find_product(rxn, rxn_flask, args)
        if rxn_flask:
            if args.verbosity > 0:
                logging.info('Products for {} have been found.'.format(rxn_condition))
            G=reaction_network(rxn_flask, tot_network,args)
            G_dict[rxn_condition]=G
        else:
            G_dict[rxn_condition]="Products are not produced."
    return G_dict

def flatten_list(lst):
    flattened = []
    for item in lst:
        if isinstance(item, list):
            flattened.extend(flatten_list(item))
        else:
            flattened.append(item)
    return flattened

def find_chemical_nodes(G):
    reaction_node, reactant_node, product_node, byproduct_node, intermediate_node, spectator_node=[], [], [], [], [], []
    mapping_number=0
    for nid, attrs in G.nodes.data():
        if nid.startswith('Reaction'):
            reaction_node.append(nid)
        elif nid.startswith('Molecule'):
#             mapping_number=max([mapping_number, max(attrs['molecule']['atom_mapping'])])
            if attrs['molecule']['identity']== 'reactant':
                reactant_node.append(nid)
            elif attrs['molecule']['identity']== 'intermediate':
                intermediate_node.append(nid)
            elif attrs['molecule']['identity']== 'product':
                product_node.append(nid)
            elif attrs['molecule']['identity']== 'byproduct':
                byproduct_node.append(nid)
            elif attrs['molecule']['identity']== 'spectator':
                spectator_node.append(nid)
    if not product_node:
        raise ValueError("Product node is not found.")
    return reaction_node, reactant_node, product_node, byproduct_node, intermediate_node, spectator_node

def create_combinations(cycle_dict):
    # Extract sequences from each group
    sequences = [seq for sequences in cycle_dict.values() for seq in sequences]
    # Create all possible combinations of one sequence from each group
    combinations = list(product(*cycle_dict.values()))
    return combinations

def sort_key(item):
    key = item[0]
    if key.startswith('Reaction'):
        return (1, int(key.split(' ')[1]))
    elif key.startswith('End'):
        return (2, 1)

def check_list_type(lst):
    if isinstance(lst, list):  # Check if it's a list
        for element in lst:
            if isinstance(element, list):  # Check for nested list
                return "Nested List"
        return "Simple List"
    else:
        return "Not a List"

def find_shared_nodes(cycles):
    shared_nodes_dict = defaultdict(list)

    if check_list_type(cycles)=="Simple List":
        shared_nodes_dict[('Single Cycle',)].append(cycles)
        return shared_nodes_dict

    # Find the common nodes for all cycle
    for i in range(len(cycles)):
        for j in range(i + 1, len(cycles)):
            common_nodes_list = list(set(cycles[i]) & set(cycles[j]))
            if common_nodes_list:
                common_nodes_list_sorted = tuple(sorted(common_nodes_list))
                if cycles[i] not in shared_nodes_dict[common_nodes_list_sorted]:
                    shared_nodes_dict[common_nodes_list_sorted].append(cycles[i])
                if cycles[j] not in shared_nodes_dict[common_nodes_list_sorted]:
                    shared_nodes_dict[common_nodes_list_sorted].append(cycles[j])

    # Check every cycle if there is a cycle without common nodes
    for cycle in cycles:
        cycle_in_shared = False
        for shared_cycles in shared_nodes_dict.values():
            if cycle in shared_cycles:
                cycle_in_shared = True
                break
        if not cycle_in_shared:
            shared_nodes_dict[('Single Cycle',)].append(cycle)

    shared_nodes_dict = {k: v for k, v in shared_nodes_dict.items() if v}

    return shared_nodes_dict

def elementary_reaction(G, args):
    '''
    Input: Reaction graph
    full: returns full reaction
    spectator: returns spectators in the reactant and product sides
    plain: returns smiles without atom mapping
    byproduct: returns byproducts in the product side
    end: returns additional reactions for the end. products>>products
    v: verbose
    '''
    # print(nx.node_link_data(G))
    if args.plain:
        smiles='smiles'
    else: smiles= 'smiles_w_mapping'
    reaction_nodes, reactant_nodes, product_nodes, byproduct_nodes, intermediate_nodes, spectator_nodes=find_chemical_nodes(G)
    elem_dict=defaultdict(lambda: list())
    if len(reaction_nodes)> args.num_reaction_node:
        raise ValueError("Too many reaction nodes.")

    # Check if there are cycles
    cycle_pathways=[i for i in nx.simple_cycles(G)]
    if args.verbosity > 2:
        logging.info('There are {} cycles'.format(len(cycle_pathways)))

    if args.verbosity > 2:
        logging.info(f'Product node is {product_nodes}')
    reaction_node_for_product=[node for node in G.predecessors(product_nodes[0])][0]
    if args.verbosity > 2:
        logging.info(f'Reaction node for product is {reaction_node_for_product}')

    if cycle_pathways:
        if len(cycle_pathways)> args.num_cycles:
            raise ValueError("Too many cycles.")
        if len(cycle_pathways)>1:
            filtered_cycles = []
            for i, cycle in enumerate(cycle_pathways):
                filtered_cycle = [node for node in cycle if node.startswith('Reaction')]
                filtered_cycles.append(filtered_cycle)
        elif len(cycle_pathways)==1:
            filtered_cycles=[node for node in cycle_pathways[0] if node.startswith('Reaction')]
        cycle_dict=find_shared_nodes(filtered_cycles)
        reaction_nodes_outside=list(set(reaction_nodes)-set(flatten_list(filtered_cycles)))
        combination_of_cycles=create_combinations(cycle_dict) # In case of more than 2 cycles
        reaction_paths=[flatten_list(list(combi)+reaction_nodes_outside) for combi in combination_of_cycles]
    else: reaction_paths = [reaction_nodes]

    # Sorting each inner list based on the numerical part of the reaction strings
    sorted_inner = [sorted(sublist, key=lambda reaction_paths: int(reaction_paths.split(' ')[1])) for sublist in reaction_paths]

    # Sorting the outer list based on the sorted inner lists
    reaction_paths = sorted(sorted_inner, key=lambda x: [int(reaction_paths.split(' ')[1]) for reaction_paths in x])

    reaction_paths = sorted(reaction_paths)
    if args.verbosity > 2:
        logging.info('There are {} paths'.format(len(reaction_paths)))

    if args.verbosity > 2:
        logging.info('Start retrieving the overall reactions')


    if args.full or args.end:
        output_smiles=list()
        for reaction_path in reaction_paths:
            # Check the reactant can be used.
            usable_reactant_nodes=list()
            produced_byproduct_nodes=list()

            for r_node in reactant_nodes:
                for reaction_node in reaction_path:
                        is_used=False
                        try:
                            paths_between_two=[path for path in nx.all_shortest_paths(G, source=r_node, target=reaction_node)]
                        except: continue
                        for path in paths_between_two:
                            path=[node for node in path if node.startswith('Reaction')]
                            not_allowed_nodes = set(path) - set(reaction_path)
                            if not not_allowed_nodes:
                                usable_reactant_nodes.append(r_node)
                                is_used=True
                            break
                        if is_used: break

            for by_node in byproduct_nodes + product_nodes:
                for reaction_node in reaction_path:
                    reaction_node_for_byproduct=[node for node in G.predecessors(by_node)][0]
                    is_produced=False
                    try:
                        paths_between_two=[path for path in nx.all_shortest_paths(G, source=reaction_node_for_byproduct, target=reaction_node)]
                    except: pass
                    for path in paths_between_two:
                        if reaction_node_for_product in path: continue
                        path=[node for node in path if node.startswith('Reaction')]
                        not_allowed_nodes = set(path) - set(reaction_path)
                        if not not_allowed_nodes:
                            produced_byproduct_nodes.append(by_node)
                            is_produced=True
                            break
                    if is_produced: break

            final_reaction_node = [node for node in G.predecessors(product_nodes[0])][0]
            if args.byproduct:
                final_products=[node for node in G.successors(final_reaction_node) if G.nodes[node]['molecule']['identity'] == 'product']
            else: final_products = product_nodes

            r_smiles=sorted([G.nodes[x]['molecule'][smiles] for x in usable_reactant_nodes])
            p_smiles=sorted([G.nodes[x]['molecule'][smiles] for x in final_products])

            not_used_reactant=list(set(reactant_nodes) - set(usable_reactant_nodes))
            if args.reagent: reagent_smiles=list()
            if args.byproduct:
                if args.reagent:
                    reagent_smiles=sorted([G.nodes[x]['molecule'][smiles] for x in produced_byproduct_nodes if x not in final_products])
                else:
                    p_smiles=p_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in produced_byproduct_nodes if x not in final_products])
            if args.spectator:
                if args.reagent:
                    reagent_smiles=reagent_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])
                else:
                    r_smiles=r_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])
                    p_smiles=p_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])

            if args.reagent:
                r_string='.'.join(r_smiles)
                reagent_string='.'.join(reagent_smiles)
                p_string='.'.join(p_smiles)
                output_smiles.append(f'{r_string}>{reagent_string}>{p_string}')
            else:
                output_smiles.append('>>'.join(['.'.join(r_smiles),'.'.join(p_smiles)]))


        if args.full: return output_smiles
        if args.end:
            for i, rxn_smi in enumerate(output_smiles):
                r_smi,reagent_smi, p_smi = rxn_smi.split('>')
                elem_dict['End'].append(f'{p_smi}>{reagent_smi}>{p_smi}')

    if args.verbosity > 2:
        logging.info('Start retrieving the elementary reactions')
    # print('reaction_paths', reaction_paths)
    for reaction_path in reaction_paths:
        for elem_node in reaction_path:
            if elem_node in elem_dict.keys(): continue
            precursor_nodes=[node for node in G.predecessors(elem_node)]
            successor_nodes=[node for node in G.successors(elem_node)]
            conmused_reactant_node=[]
            # print('elem_node', elem_node, precursor_nodes)
            if args.spectator:
                for r_node in reactant_nodes:
                    for precursor in precursor_nodes:
                        is_used=False
                        if precursor in reactant_nodes:
                            conmused_reactant_node.append(precursor)
                            is_used=True
                            continue
                        try:
                            paths_between_two=[path for path in nx.all_shortest_paths(G, source=r_node, target=precursor)]
                        except: continue
                        for path in paths_between_two:
                            path=[node for node in path if node.startswith('Reaction')]
                            not_allowed_nodes = set(path) - set(reaction_path)
                            if not not_allowed_nodes:
                                conmused_reactant_node.append(r_node)
                                is_used=True
                                break
                        if is_used: break

                not_used_reactant=list(set(reactant_nodes) - set(conmused_reactant_node))
            else: not_used_reactant=[]

            if args.byproduct:
                produced_byproduct_nodes=list()
                reaction_node_for_product = [node for node in G.predecessors(product_nodes[0])][0]
                for by_node in byproduct_nodes + product_nodes:
                    # print('by_node is ', by_node)
                    reaction_node_for_byproduct=[node for node in G.predecessors(by_node)][0]
                    if reaction_node_for_byproduct not in reaction_path:
                        # print(f'{reaction_node_for_byproduct} node is not in the path')
                        # print(f'reaction_path is {reaction_path}')
                        continue
                    # print('0000',by_node, reaction_node_for_byproduct)
                    for precursor in precursor_nodes:
                        # print('precursor', precursor)
                        is_produced=False
                        if precursor in reactant_nodes:
                            # print('This is reactant node!')
                            continue
                        try:
                            paths_between_two=[path for path in nx.all_shortest_paths(G, source=reaction_node_for_byproduct, target=precursor)]
                        except Exception as e:
                            # print(e)
                            continue
                        # print('hi1', paths_between_two)

                        for path in paths_between_two:
                            # print('reaction_node_for_product', reaction_node_for_product)
                            # print('path', path)
                            if reaction_node_for_product not in path:
                                # print('hi2')
                                is_there_reactant = any(element in reactant_nodes for element in path)
                                if is_there_reactant: continue
                                if reaction_node_for_product in path:
                                    further_reacting_intermidiate = [node for node in
                                                                     G.successors(reaction_node_for_product) if
                                                                     G.nodes[node]['molecule'][
                                                                         'identity'] == 'intermediate']
                                    if not further_reacting_intermidiate:
                                        continue
                                path=[node for node in path if node.startswith('Reaction')]
                                not_allowed_nodes = set(path) - set(reaction_path)
                                if not not_allowed_nodes and elem_node not in path:
                                    produced_byproduct_nodes.append(by_node)
                                    is_produced=True
                                    break
                        if is_produced: break
            else: produced_byproduct_nodes=list()
            r_smiles=sorted([G.nodes[x]['molecule'][smiles] for x in precursor_nodes])
            p_smiles=sorted([G.nodes[x]['molecule'][smiles] for x in successor_nodes])

            if args.reagent:
                reagent_smiles=list()

            if args.byproduct:
                if args.reagent:
                    reagent_smiles=reagent_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in produced_byproduct_nodes])
                else:
                    r_smiles=r_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in produced_byproduct_nodes])
                    p_smiles=p_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in produced_byproduct_nodes])
            if args.spectator:
                if args.reagent:
                    reagent_smiles=reagent_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])
                else:
                    r_smiles=r_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])
                    p_smiles=p_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])

            if args.reagent:
                r_string='.'.join(r_smiles)
                reagent_string='.'.join(reagent_smiles)
                p_string='.'.join(p_smiles)
                elem_dict[elem_node].append(f'{r_string}>{reagent_string}>{p_string}')
            else:
                elem_dict[elem_node].append('>>'.join(['.'.join(r_smiles),'.'.join(p_smiles)]))

    elem_dict = dict(sorted(elem_dict.items(), key=sort_key))
    elementary_reaction_output=list()

    for elem_rxn_smi in flatten_list([value for value in elem_dict.values()]):
        if elem_rxn_smi not in elementary_reaction_output:
            elementary_reaction_output.append(elem_rxn_smi)
    return elementary_reaction_output
