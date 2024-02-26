import os
import sys
import json
import rdkit
import gzip
from rdkit import Chem
from rdkit.Chem import AllChem
from tqdm import tqdm
from joblib import Parallel, delayed
import itertools 
from itertools import product, permutations
from templates import acidBase_lookup, reaction_templates
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

def collect_reaction_class(reactions):
    '''
    To collect one reaction class from list of reaction, in the form of {'rxnsmiles': {'reaction_name': name}}.
    reactions: list of reactions 
    lst_rxn: A reaction list
    rxn_class_name: Reaction name to be collected
    '''
    rxns_grouped_by_class = {}
    unclassified = []
    
    for reaction in reactions:
        if 'reaction_name' in reaction.keys():
            rxns_grouped_by_class[reaction['reaction_name']] = rxns_grouped_by_class.get(reaction['reaction_name'], []) + [reaction]
        else: 
            print("No reaction name for", reaction)
            unclassified.append(reaction)
        
    
    print(len(unclassified))
    
    return rxns_grouped_by_class

def get_class_key(class_name):
    
    '''
    returns the key(tuple) of class_reaction_templates that includes the class_name
    '''
    
    for classes in Reaction_templates.class_reaction_templates.keys():
        if class_name in classes:
            return classes
        
    return None
 

def reagent_matching_for_class(class_name, reactions):
    

    matched_reactions = {}
    class_key = get_class_key(class_name)
    
    if class_key:
        results = Parallel(n_jobs=10)(delayed(reagent_matching_for_single_reaction)(class_key, reaction) for reaction in tqdm(reactions))
        return results
    
    else:
        print("No class key for class", class_name)
        return None

    
def reagent_matching_for_single_reaction(class_key, reaction):
    '''
    To match reagents for a certain class from list of reaction, in the form of {'rxnsmiles': {'reaction_name': name}}.
    reactions: list of reactions 
    '''
    reactants, agents, products = reaction['reaction_smiles'].split(">")
    mols = [Chem.MolFromSmiles(smi) for smi in reactants.split(".")+agents.split(".")]
    
    class_data = Reaction_templates.class_reaction_templates[class_key]
    
    reaction['conditions'] = []
    
    for cond_name, cond_data in class_data.items():
        # print(cond_name, cond_data)
        if cond_data['Reagent']:
            cond_mols = [Chem.MolFromSmarts(x) for x in cond_data['Reagent']]
            matched_reagents = []
            for patt in cond_mols:
                for mol in mols:
                    if mol.GetSubstructMatch(patt):
                        matched_reagents.append(mol)
                        break

            if len(cond_mols) == len(matched_reagents):     
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

def prepare_reactants(reaction_dict):
    """
    Preprocessing to make reactants with isotope labeling
    reaction_dict: one reaction in reactions_with_conditions    
    """
    
    rxn_smi=reaction_dict['reaction_smiles']

    reactants, agents, products = [mols_from_smiles_list(x) for x in 
                                [mols.split('.') for mols in rxn_smi.split('>')]]
    
    rmol=reactants+agents
    
    idx=1
    reactant_pool=list()
    reactant_dict=dict()
    _smiles_list=[]
    for num, mol in enumerate(rmol):
        start_idx=idx
        for atom in mol.GetAtoms():
            atom.SetIsotope(idx) 
            idx+=1
            
        _smiles_list.append(Chem.MolToSmiles(mol))
        
        #Save molecule in a dictionary
        
        reactant_dict[num] = {'smiles_w_isotope':Chem.MolToSmiles(mol),   # Using isotope as atom mapping for reaction
                            'atom_mapping':[i for i in range(start_idx,idx)],  #Checking for atom-mapping collision
                            'smiles':Chem.MolToSmiles(remove_atom_map(mol), isomericSmiles=False),  #Plain SMILES string
                           'smiles_w_mapping': Chem.MolToSmiles(isotope_to_atommap(mol), isomericSmiles=False), # Output for generating elementary reaction graph
                           'rxn_history': list(),
                            'identity': 'reactant'
                             }
        
    reactant_dict['last_mapping_number']=idx-1
    return reactant_dict

def check_products_validity(outcome):
 #Many different outcomes can be producted from one template (e.g. symmetric molecules)

    good_molecule=0 # To check every molecules make sense
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

def find_reactants(rxn_flask, rxn_templates, stoichiometry):
    reactive_dict = dict()
    for templ in rxn_templates:
        reactant_dict=dict()
        r, p = templ.split('>>')
#         if '.' in r:
        r=r.split('.')
#         else: 
#             r=[r]        
        num_reactants=len(r)    

        mol_dict = {reactant: Chem.MolFromSmiles(rxn_flask[reactant]['smiles']) for reactant in rxn_flask}
        templ_dict = {reac_temp: Chem.MolFromSmarts(reac_temp) for reac_temp in r}

        # Find reactive chemical species
        for mol_key in mol_dict:        
            mol = mol_dict[mol_key]
            for temp_key in templ_dict:
                pat = templ_dict[temp_key]
                if mol.GetSubstructMatch(pat):
                    if temp_key not in reactant_dict:
                        reactant_dict[temp_key] = []
                    reactant_dict[temp_key].append(mol_key)
        print(reactant_dict)
        combinations = [list(combination) for combination in itertools.product(*reactant_dict.values()) if len(combination)==num_reactants]
        
        # Check the atom map collision
        reactive_combination = []
        for combination in combinations:
            reactant_history_per_com=list()
            if num_reactants>1: 
                reactant_atommap=[set(i) for i in [rxn_flask[num]['atom_mapping'] for num in combination]]
                # Check if there is the same atom in two reactants
                if reactant_atommap[0]&reactant_atommap[1]: 
                    continue
                reactant_history=[i for i in [rxn_flask[num]['rxn_history'] for num in combination]]
                reactant_history.extend(combination)
                
                # Check if two species come from the same reactant
                if has_duplicates(reactant_history): 
                    continue        
            reactive_combination.append(combination)
            reactant_history_per_com.append(reactant_history)
            if num_reactants>1:
                reactive_combination.append(list(reversed(combination)))
                reactant_history_per_com.append(reactant_history)
            reactive_dict[templ]={'combination': reactive_combination,
                                  'num_reactants': num_reactants,
                                  'reactant_history': reactant_history_per_com}
            
    return reactive_dict

def remove_isotope(prod_mol):
    prod_istope_smi=Chem.MolToSmiles(prod_mol)
    for a in prod_mol.GetAtoms():
        a.SetIsotope(0)
    return Chem.MolFromSmiles(prod_istope_smi), Chem.MolToSmiles(prod_mol)

def product_dict(prod_mol, reactant_id, reactant_history):
    prod_smi_isotope = Chem.MolToSmiles(prod_mol) # Get isotompe-mapped smiles
    atom_mapping=list()

    for a in prod_mol.GetAtoms(): #Change isotope to atom map
        a.SetAtomMapNum(a.GetIsotope())
        atom_mapping.append(a.GetIsotope())
        a.SetIsotope(0)
    prod_smi_map=Chem.MolToSmiles(prod_mol) #Get atom mapped smiles
    for a in prod_mol.GetAtoms(): #Remove any mapping
        a.SetAtomMapNum(0)
    prod_smi=Chem.MolToSmiles(prod_mol) #Get plain smiles

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

def proton_balanced_template(rxn_flask, pKas, rxn_templates):
    
    """
    rxn_flask: dictionary - key: species index (0~), value: dictionary of species data (smiles, tags, etc)
    pKas: list of pKas for each rxn_template in rxn_templates
    """
    proton_balanced_rxn_templates=list()
    if pKas == [None]*len(pKas):
#         for data in rxn_flask.values():
#             denotes whether a species works as acid or base in each template
#             data['acid_or_base'] = [None for temp in rxn_templates] 
#         return rxn_flask, rxn_templates
        return rxn_templates
    
    for pKa in pKas:
        if not pKa: continue
        A=pKa.get('A')
        B=pKa.get('B')      
        if A:
            filtered_data = [d for d in AcidBase_lookup.Acid_base if 'A' in d['role'] and d['pKa'] < A]
            sorted_data = sorted(filtered_data, key=lambda x: x['pKa'])
            possible_acid_base=find_acid_base(rxn_flask, sorted_data, 'A')
        if B:
            filtered_data = [d for d in AcidBase_lookup.Acid_base if 'B' in d['role'] and d['pKa'] > B]
            sorted_data = sorted(filtered_data, key=lambda x: x['pKa'], reverse=True)
            possible_acid_base=find_acid_base(rxn_flask, sorted_data, 'B')
        if not possible_acid_base: continue
        

        for acid_base in possible_acid_base:
            for templ in rxn_templates:
                r, p = templ.split('>>')
                new_r = '.'.join([r, acid_base[0]])
                new_p = '.'.join([p, acid_base[1]])
                templ='>>'.join([new_r,new_p])
                proton_balanced_rxn_templates.append(templ)       
        
    # return rxn_flask, proton_balanced_rxn_templates
    return proton_balanced_rxn_templates

def allow_unimolecular_rxn(rxn_templates):
    new_templates=[]
    for templ in rxn_templates:
        rxn = AllChem.ReactionFromSmarts(templ)
        num_reactants=rxn.GetNumReactantTemplates()
        
        new_templates.append(templ)
        if num_reactants > 1:
            r, p = templ.split('>>')
            r= '('+r+')'
            p='('+p+')'
            new_temp='>>'.join([r,p])
            new_templates.append(new_temp)
    return new_templates
    

def find_acid_base(rxn_flask, filtered_list, ab_condition):
    possible_acid_base=[]
    for acid_base in filtered_list:
        if ab_condition=='A':
            reactant=acid_base['Acid']
            product=acid_base['Base']        
        elif ab_condition=='B':
            
            reactant=acid_base['Base']
            product=acid_base['Acid']        
        
        patt = Chem.MolFromSmarts(reactant)
        mols = [Chem.MolFromSmiles(rxn_flask[x]['smiles']) for x in rxn_flask  if type(x) is int]
        for mol in mols:
            if mol and mol.GetSubstructMatch(patt):
                possible_acid_base.append([reactant, product])
    return possible_acid_base

def find_reactants(rxn_flask, rxn_templates, stoichiometry, v=False):
    reactive_dict = dict()
    mol_dict = {reactant: Chem.MolFromSmiles(rxn_flask[reactant]['smiles']) for reactant in rxn_flask if type(reactant) is int}
    
    for num_temp, templ in enumerate(rxn_templates):
        reactant_dict=dict()
        
        rxn = AllChem.ReactionFromSmarts(templ)
        r = [rmol for rmol in rxn.GetReactants()]
        num_reactants=rxn.GetNumReactantTemplates()
        
        templ_dict = {templ: r}
        # Find reactive chemical species
        for mol_key in mol_dict:        
            mol = mol_dict[mol_key]
            for temp_key in templ_dict:
                patterns = templ_dict[temp_key]
                for pat in patterns:
                    if mol and mol.GetSubstructMatch(pat):
                        if Chem.MolToSmarts(pat) not in reactant_dict:
                            reactant_dict[Chem.MolToSmarts(pat)] = []
                        reactant_dict[Chem.MolToSmarts(pat)].append(mol_key)
        
        if v and reactant_dict: 
            print('A dictionary of all possible reactants for a template of {}'.format(templ))
            print(reactant_dict, '\n')
            
        combinations = [list(combination) for combination in itertools.product(*reactant_dict.values()) if len(combination)==num_reactants]
        # Check the atom map collision
        reactive_combination = list()
        reactant_history_per_com=list()
        if stoichiometry:
            # To save every combination
            duplicated_combi=list()
            
        for combination in combinations:            
            if num_reactants>1: 
                reactant_atommap=[set(i) for i in [rxn_flask[mol_num]['atom_mapping'] for mol_num in combination]]
                # Check if there is the same atom in two reactants
                if reactant_atommap[0]&reactant_atommap[1]: 
                    continue
            reactant_history=[i for i in [rxn_flask[mol_num]['rxn_history'] for mol_num in combination]]
            reactant_history.extend(combination)
            # Check if at least one species came from the same reactant
            if has_duplicates(reactant_history): 
                if stoichiometry:
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
            
        if stoichiometry and duplicated_combi and not reactive_dict:
            idx = rxn_flask['last_mapping_number']+1
            new_combination_list=list()
            
            for combination in duplicated_combi:
                new_combination=list()
                for mol_num in combination:
                    if rxn_flask[mol_num]['identity']=='reactant':
                        start_idx=idx
                        new_mol = Chem.MolFromSmiles(rxn_flask[mol_num]['smiles'])
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
                        new_combination.append(len(rxn_flask))
                    else: new_combination.append(mol_num)
                new_combination_list.append(new_combination)
                        
            rxn_flask['last_mapping_number']=idx-1
            
            for combination in new_combination_list:
                reactant_history_per_com=list()
                if num_reactants>1: 
                    reactant_atommap=[set(i) for i in [rxn_flask[mol_num]['atom_mapping'] for mol_num in combination]]
                    # Check if there is the same atom in two reactants
                    if reactant_atommap[0]&reactant_atommap[1]: 
                        continue
                reactant_history=[i for i in [rxn_flask[mol_num]['rxn_history'] for mol_num in combination]]
                reactant_history.extend(combination)
                # Check if at least one species came from the same reactant
                if has_duplicates(reactant_history):
                    continue        

                reactive_combination.append(combination)
                reactant_history_per_com.append(reactant_history)

                if num_reactants>1:
                    reactive_combination.append(list(reversed(combination)))
                    reactant_history_per_com.append(reactant_history)
                reactive_dict[templ]={'combination': reactive_combination,
                                      'num_reactants': num_reactants,
                                      'reactant_history': reactant_history_per_com}
            
    if v and not reactive_dict:
        print('No reactants are found!')
    return rxn_flask, reactive_dict

def has_duplicates(input_list):
    lst=list()
    elem=list()
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



def run_single_reaction(rxn_flask, single_step, proton=False, uni_rxn=False, stoichiometry=False, v=False):
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
    
    reaction_network=[]
    reaction_pair=[]
    
    if v: 
        num_original_temp=len(rxn_templates)
    
    if proton: # Get proton balanced reaction template
        # TODO: label acids and bases (e.g. 'acid' instead of just 'reactant' when labelling), add isAcid, isBase
        rxn_templates=proton_balanced_template(rxn_flask, pKas, rxn_templates)
        num_proton_temp=len(rxn_templates)
        if v:
            print('New proton balanced {} templates are generated'.format(num_proton_temp-num_original_temp))
            num_original_temp=len(rxn_templates)
        
    if uni_rxn: # Get unimolecular reaction template
        rxn_templates=allow_unimolecular_rxn(rxn_templates)
        if v:
            num_uni_temp=len(rxn_templates)
            print('New unimolecular {} templates are generated \n'.format(num_uni_temp-num_original_temp))
    
    #Define the reactive chemicals
    #rxn_flask may be updated if stoichiometry is on
    rxn_flask, reactive_dict=find_reactants(rxn_flask, rxn_templates, stoichiometry, v)
    
    for templ in reactive_dict:        
        if v: print('Template is ', templ)
        
        combinations=reactive_dict[templ]['combination']
        num_reactants=reactive_dict[templ]['num_reactants']
        reactant_history=reactive_dict[templ]['reactant_history']
        
        rxn = AllChem.ReactionFromSmarts(templ)
        for combination, history in zip(combinations, reactant_history):                             
            reactants=[Chem.MolFromSmiles(rxn_flask[num]['smiles_w_isotope'],sanitize=False)
                                 for num in combination]  
            outcomes = rxn.RunReactants(reactants)
            if not outcomes:
                continue
                
            if v:
                reactants_smi=[rxn_flask[num]['smiles'] for num in combination]  
                print('The reactants are {}'.format(reactants_smi))
                
            for outcome in outcomes:
                reactant_id=set(num for num in combination)
                if not check_products_validity(outcome): continue
                product_id=set()
                if v: print('{} products are formed \n'.format(len(outcome)))
                for prod_mol in outcome:     
                    new_prod=product_dict(prod_mol, reactant_id, history)
                    r_map = set(flatten_list([rxn_flask[rid]['atom_mapping'] for rid in reactant_id]))
                    prod_num = [key for key, value in rxn_flask.items() if type(key) is int and value['smiles'] == new_prod['smiles'] and sorted(value['atom_mapping']) == sorted(new_prod['atom_mapping'])]
                    p_map = set(flatten_list([rxn_flask[rid]['atom_mapping'] for rid in prod_num]))
                    
                    if prod_num and r_map.issuperset(p_map): 
                        product_id.add(prod_num[0])
                        continue
                    elif prod_num and not r_map.issuperset(p_map): 
                        continue
                    else: 
                        num_mols=len(rxn_flask)+1
                        product_id.add(num_mols)
                        rxn_flask[num_mols]= new_prod
                        
                if [reactant_id,product_id] not in reaction_pair:
                    reaction_pair.append([reactant_id,product_id])
                    reaction_network.append([reactant_id,templ, product_id])
                
    return rxn_flask, reaction_network

def run_full_reaction(rxn_flask, condition, v=False):
    
    """
    Given an initial rxn_flask, applies the templates for each condition and
    returns lists of final reaction flasks and networks for each condition
    Returns empty lists if there are no conditions
    
    rxn_flasks: a list of final reaction flasks for each condition
    tot_networks: a list of tot_net for each condition
    """
    
    
    tot_networks = list()
    rxn_flasks = list()
    for cond in condition:
        tot_net = list()
        rxn_flask_cond = rxn_flask
        for temp in cond['Stages'].values():
            rxn_flask_cond, reaction_network=run_single_reaction(rxn_flask_cond , temp, uni_rxn=True, proton=True, stoichiometry=True, v=v)        
            tot_net+=reaction_network
            
        tot_networks.append(tot_net)
        rxn_flasks.append(rxn_flask_cond)
            
            
    return rxn_flasks, tot_networks

def find_product(example_rxn, rxn_flask):
    reaction_smi = example_rxn['reaction_smiles']
    react_smi, reagent_smi, product_smi = reaction_smi.split('>')
    
    product_smi_list=product_smi.split('.')
    real_product_smi_list=[]
    for psmi in product_smi_list:
        pmol=Chem.MolFromSmiles(psmi)
        pat = Chem.MolFromSmarts("[#6]")
        if len(pmol.GetSubstructMatches(pat)) > 0:
            pmol=remove_atom_map(pmol)
            real_product_smi_list.append(Chem.MolToSmiles(pmol))
    
    matching_keys = [key for key, value in rxn_flask.items() if type(key) is int and value['smiles'] in real_product_smi_list]
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

def reaction_network(rxn_flask, tot_network, simple=False, light=True, full=False):
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
            if light:
                molecule={'identity': rxn_flask[chemical]['identity'],
                         'smiles_w_mapping': rxn_flask[chemical]['smiles_w_mapping'],
                         'smiles': rxn_flask[chemical]['smiles']}
            else: molecule=rxn_flask[chemical]
            G.add_node(chemical_node, molecule=molecule)
            G.add_edge(chemical_node, reaction_node)            
        for chemical in last_set:
            chemical_node=f"Molecule {chemical}"
            if light:
                molecule={'identity': rxn_flask[chemical]['identity'],
                         'smiles_w_mapping': rxn_flask[chemical]['smiles_w_mapping'],
                         'smiles': rxn_flask[chemical]['smiles']}
            else: molecule=rxn_flask[chemical]
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
            
    # Find all the routes connecting the reactants and products
    reaction_path=list()
    for rid in reactant_ids:
        for pid in product_ids:
            if simple:
                if nx.all_simple_paths(G, source=rid, target=pid):
                    for path in nx.all_simple_paths(G, source=rid, target=pid):
                        if len(reaction_nodes_in_path) < 7:
                            reaction_path.append(node for node in path if node.startswith('Reaction'))
            else:
                if nx.all_shortest_paths(G, source=rid, target=pid):
                    for path in nx.all_shortest_paths(G, source=rid, target=pid):
                        reaction_nodes_in_path=[node for node in path if node.startswith('Reaction')]
                        if len(reaction_nodes_in_path) < 7:
                            reaction_path.append([node for node in path if node.startswith('Reaction')])

    
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
#     print(path_node)
    if full: return G
    
    # Make subgraph of a true path    
    G_sub = nx.DiGraph(G.subgraph(path_node))
    
    #Check if there is not-connected intermediates.
    not_connected_inter_nodes=list()
    for nid, attrs in G_sub.nodes.data():
        if nid.startswith('Molecule') and attrs['molecule']['identity']== 'intermediate':
            parents=[node for node in G_sub.predecessors(nid)]
            child=[node for node in G_sub.successors(nid)]
            if not parents or not child:
                not_connected_inter_nodes.append(nid)
    if not_connected_inter_nodes:
        for rid, pid in permutations(not_connected_inter_nodes, 2):
            try:
                shortest_paths = list(nx.all_shortest_paths(G, source=rid, target=pid))
                if shortest_paths:
                    for path in shortest_paths:
                        reaction_path.append([node for node in path if node.startswith('Reaction')])
                    reaction_path=flatten_list(reaction_path)
                    for reaction_node in reaction_path:
                        path_node.append(reaction_node)
                        neighbor_node=[n for n in nx.all_neighbors(G, reaction_node)]
                        for nn in neighbor_node:
                            path_node.append(nn)
                    path_node=list(set(path_node))
                        
                    G_sub = nx.DiGraph(G.subgraph(path_node))
            except nx.NetworkXNoPath: continue
                
    
    # Add spectators
    missing_molecule_nodes = [node_id for node_id in reactant_list if f'Molecule {node_id}' not in G.nodes]
    
    for node_id in missing_molecule_nodes:
        chemical_node=f"Molecule {node_id}"
        if light:
            molecule={'identity': rxn_flask[node_id]['identity'],
                     'smiles_w_mapping': rxn_flask[node_id]['smiles_w_mapping'],
                     'smiles': rxn_flask[node_id]['smiles']}
        else: molecule=rxn_flask[node_id]
                
        molecule['identity']='spectator'
        G_sub.add_node(chemical_node, molecule=molecule)
    
    return G_sub

def draw_reaction_graph(G, size=[5,5], labels = True):
#     G=nx.node_link_graph(G)
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

def get_mechanistic_network(rxn, v=False, simple=False, light=True):
    condition = calling_rxn_template(rxn)
    rxn_flask=prepare_reactants(rxn)
    rxn_flasks, tot_networks=run_full_reaction(rxn_flask, condition, v=v)
    G_list=list()
    for rxn_flask,tot_network in zip(rxn_flasks, tot_networks):
        rxn_flask=find_product(rxn, rxn_flask)
        if rxn_flask:
            G=reaction_network(rxn_flask, tot_network, simple=simple, light=light)
            G_list.append(G)
            
    if G_list:
        return G_list
    else: 
        raise ValueError("Products are not produced.")
        
        
        
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
    
    # Check every cycles if there is a cycle without common nodes
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

def elementary_reaction(G, full=False, spectator=True, plain=False, byproduct=True, end=True, reagent=False, v=False):
    '''
    Input: Reaction graph
    full: returns full reactions
    spectator: returns spectators in the reactant and product sides
    plain: returns smiles without atom mapping
    byproduct: returns byproducts in the product side
    end: returns additional reactions for the end. products>>products
    v: verbose
    '''
    if plain:
        smiles='smiles'
    else: smiles= 'smiles_w_mapping'
   
    reaction_nodes, reactant_nodes, product_nodes, byproduct_nodes, intermediate_nodes, spectator_nodes=find_chemical_nodes(G)
#     elem_dict=dict()
    elem_dict=defaultdict(lambda: list())
    # Check if there are cycles
    cycle_pathways=[i for i in nx.simple_cycles(G)]
    
    reaction_node_for_product=[node for node in G.predecessors(product_nodes[0])][0]
    
    if cycle_pathways:
        if len(cycle_pathways)>1: 
            filtered_cycles = []
            for cycle in cycle_pathways:
                filtered_cycle = [node for node in cycle if node.startswith('Reaction')]
                filtered_cycles.append(filtered_cycle)
        elif len(cycle_pathways)==1: 
            filtered_cycles=[node for node in cycle_pathways[0] if node.startswith('Reaction')]
        cycle_dict=find_shared_nodes(filtered_cycles)
        reaction_nodes_inside=list(set(flatten_list(filtered_cycles)))
        
        reaction_nodes_outside=list(set(reaction_nodes)-set(flatten_list(filtered_cycles)))
        combination_of_cycles=create_combinations(cycle_dict) # In case of more than 2 cycles
        reaction_paths=[flatten_list(list(combi)+reaction_nodes_outside) for combi in combination_of_cycles]
    else: reaction_paths = [reaction_nodes]
        
    if full or end:
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
                        
            for by_node in byproduct_nodes:
                for reaction_node in reaction_path:
                    reaction_node_for_byproduct=[node for node in G.predecessors(by_node)][0]
                    is_produced=False
                    try:
                        paths_between_two=[path for path in nx.all_shortest_paths(G, source=reaction_node_for_byproduct, target=reaction_node)]
                    except: pass
                    for path in paths_between_two:
                        if reaction_node_for_product in path: continue
                        path=[node for node in path if node.startswith('Reaction')]

                        if not not_allowed_nodes:
                            produced_byproduct_nodes.append(by_node)
                            is_produced=True
                            break
                    if is_produced: break
                        
            final_reaction_node = [node for node in G.predecessors(product_nodes[0])][0]
            if byproduct:
                final_products=[node for node in G.successors(final_reaction_node)]
            else: final_products = product_nodes
            
            r_smiles=sorted([G.nodes[x]['molecule'][smiles] for x in usable_reactant_nodes])
            p_smiles=sorted([G.nodes[x]['molecule'][smiles] for x in final_products])
            
            not_used_reactant=list(set(reactant_nodes) - set(usable_reactant_nodes))
            if reagent: reagent_smiles=list()
            if byproduct:
                if reagent:
                    reagent_smiles=sorted([G.nodes[x]['molecule'][smiles] for x in produced_byproduct_nodes if x not in final_products])
                else:
                    p_smiles=p_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in produced_byproduct_nodes if x not in final_products])
            if spectator:
                if reagent:
                    reagent_smiles=reagent_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])
                else:
                    r_smiles=r_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])
                    p_smiles=p_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])
            
            if reagent:
                r_string='.'.join(r_smiles)
                reagent_string='.'.join(reagent_smiles)
                p_string='.'.join(p_smiles)
                output_smiles.append(f'{r_string}>{reagent_string}>{p_string}')
            else:
                output_smiles.append('>>'.join(['.'.join(r_smiles),'.'.join(p_smiles)]))
 
                
        if full: return output_smiles
        if end: 
            for i, rxn_smi in enumerate(output_smiles):
                r_smi,reagent_smi, p_smi = rxn_smi.split('>')
                elem_dict['End'].append(f'{p_smi}>{reagent_smi}>{p_smi}')
                
    for reaction_path in reaction_paths:
        for elem_node in sorted(reaction_path):
            precursor_nodes=[node for node in G.predecessors(elem_node)]
            successor_nodes=[node for node in G.successors(elem_node)]
            conmused_reactant_node=[]

            if spectator: 
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

            if byproduct:  
                produced_byproduct_nodes=list()
                for by_node in byproduct_nodes:
                    reaction_node_for_byproduct=[node for node in G.predecessors(by_node)][0]
                    for precursor in precursor_nodes:
                        is_produced=False
                        if precursor in reactant_nodes:
                            continue                            
                        try:
                            paths_between_two=[path for path in nx.all_shortest_paths(G, source=reaction_node_for_byproduct, target=precursor)]
                            for path in paths_between_two:

                                if reaction_node_for_product in path: continue
                                is_there_reactant = any(element in reactant_nodes for element in path)
                                if is_there_reactant: continue
                                path=[node for node in path if node.startswith('Reaction')]
                                not_allowed_nodes = set(path) - set(reaction_path)
                                if not not_allowed_nodes and elem_node not in path:
                                    produced_byproduct_nodes.append(by_node)
                                    is_produced=True
                                    break
                            if is_produced: break
                        except Exception as e:
                            if v:
                                print(e)
                            pass
            else: produced_byproduct_nodes=list()

            r_smiles=sorted([G.nodes[x]['molecule'][smiles] for x in precursor_nodes])
            p_smiles=sorted([G.nodes[x]['molecule'][smiles] for x in successor_nodes])

            if reagent:
                reagent_smiles=list()

            if byproduct:
                if reagent: 
                    reagent_smiles=reagent_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in produced_byproduct_nodes])
                else:
                    r_smiles=r_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in produced_byproduct_nodes])
                    p_smiles=p_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in produced_byproduct_nodes])
            if spectator:
                if reagent:
                    reagent_smiles=reagent_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])
                else:
                    r_smiles=r_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])
                    p_smiles=p_smiles+sorted([G.nodes[x]['molecule'][smiles] for x in spectator_nodes+not_used_reactant])

            if reagent:
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
    