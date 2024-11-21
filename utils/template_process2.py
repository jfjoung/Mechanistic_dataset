from templates import AcidBase_lookup
import itertools
import logging
import networkx as nx
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from utils.exceptions import *
from utils.explicit_H import  modify_explicit_H
RDLogger.DisableLog('rdApp.*')

def flatten_list(lst):
    flattened = []
    seen = set()  # To keep track of unique elements

    for item in lst:
        if isinstance(item, list):
            # Recursively flatten and deduplicate
            for sub_item in flatten_list(item):
                if sub_item not in seen:
                    flattened.append(sub_item)
                    seen.add(sub_item)
        else:
            if item not in seen:
                flattened.append(item)
                seen.add(item)

    return flattened
def remove_atom_map(mol, isotope=False):
    '''
    Removing the atom mapping and isotope labeling if necessary
    '''
    mol_copy = Chem.Mol(mol)
    for atom in mol_copy.GetAtoms():
        atom.SetAtomMapNum(0)
        if isotope: atom.SetIsotope(0)
    return mol_copy

def isotope_to_atommap(mol, get_map_list = False):
    mol_copy = Chem.Mol(mol)
    map_list = []
    for idx, a in enumerate(mol_copy.GetAtoms()):
        if get_map_list:
            map_list.append(a.GetIsotope())
        a.SetAtomMapNum(a.GetIsotope())
        a.SetIsotope(0)
    if get_map_list:
        return mol_copy, map_list
    else: return mol_copy

def get_plain_smiles(mol):
    '''
    Convert mol object to plain SMILES, no atom mapping and no isotope labeling
    '''
    _mol = Chem.Mol(mol)
    plain = Chem.MolToSmiles(remove_atom_map(_mol, isotope=True))# , isomericSmiles=False)
    return plain
    
def has_duplicate_atom(atom_map_list):
    new_map_list = flatten_list(atom_map_list)
    if len(new_map_list) != len(set(new_map_list)):
            return True
    return False

def get_max_atom_map(mol):
    max_mapping = max(atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0)
    return max_mapping


class Template_process:
    '''
    To process anything related to templates.
    Before enumerating the templates, find if the reactants are present in flask.
    '''

    def __init__(self, template_dict, args):
        self.template_list = template_dict['Templates']
        self.pKa_list = template_dict['pKa']
        self.description = template_dict['Description']
        self.args = args

        if self.args.explicit_H:
            self.ps = Chem.SmilesParserParams()
            self.ps.removeHs = False
            self.ps.sanitize = False
        else:
            self.ps = Chem.SmilesParserParams()
            self.ps.sanitize = False
        
        if self.args.verbosity:
            logging.info(f'Templates of elementary step for {self.description}')

    def __str__(self):
        return f'Templates of elementary step for {self.description}'

    def check_num_templates(self):
        if len(self.template_list) > self.args.max_num_temp:
            raise TooManyTemplates

    def find_reactants(self, node, reactant_node):

        self.check_num_templates()

        template_reactant_dict = {} #Get every reactant for each template
        for idx, templ in enumerate(self.template_list):
            # print(templ)
            rxn = AllChem.ReactionFromSmarts(templ)
            patterns = [rmol for rmol in rxn.GetReactants()] #Get reactants as described in the templates

            templ_mol_pair = {}  #Save the pattern - molecule pair
            num_match = 0

            for pat in patterns:
                pat.UpdatePropertyCache(strict=False)
                # print('Pat :', Chem.MolToSmarts(pat))
                possible_reactant_list = []  #Save the matched reactant
                for mol in node['mol']:
                    mol.UpdatePropertyCache(strict=False)
                    # print(Chem.MolToSmiles(mol))
                    if mol and mol.GetSubstructMatch(pat): # and mol_node.smiles not in possible_reactant_smiles_list:
                        possible_reactant_list.append(mol)
                        # print('Mol :', Chem.MolToSmiles(mol), 'Pat :', Chem.MolToSmarts(pat))

                if possible_reactant_list:
                    num_match += 1
                    templ_mol_pair[pat] = possible_reactant_list

            if num_match == len(patterns):
                template_reactant_dict[templ] = templ_mol_pair

        # print(template_reactant_dict)

        if self.args.stoichiometry and template_reactant_dict == {}:
            # print(num_match, 'templ', templ)
            # print(reactant_node)
            # print('templ_mol_pair', [pat for pat in templ_mol_pair.keys()])
            max_atom_map = get_max_atom_map(Chem.MolFromSmiles(node['smiles_w_mapping'], sanitize=False)) + 1
            for idx, templ in enumerate(self.template_list):
                rxn = AllChem.ReactionFromSmarts(templ)
                patterns = [rmol for rmol in rxn.GetReactants()]
                templ_mol_pair = {}
                num_match = 0
                for pat in patterns: 
                    possible_reactant_list = []

                    for mol in node['mol']:
                        mol.UpdatePropertyCache(strict=False)
                    # print(Chem.MolToSmiles(mol))
                    if mol and mol.GetSubstructMatch(pat): # and mol_node.smiles not in possible_reactant_smiles_list:
                        possible_reactant_list.append(mol)
                        # print('Mol :', Chem.MolToSmiles(mol), 'Pat :', Chem.MolToSmarts(pat))

                        if possible_reactant_list:
                            num_match += 1
                            templ_mol_pair[pat] = possible_reactant_list
                    if pat not in templ_mol_pair.keys():
                        # print('Pat :', Chem.MolToSmarts(pat))
                        for mol in reactant_node['mol']:
                            if mol and mol.GetSubstructMatch(pat):
                                # print('After adding Mol :', Chem.MolToSmiles(mol), 'Pat :', Chem.MolToSmarts(pat))
                                new_smi = Chem.MolToSmiles(mol)
                                # print('new_smi', new_smi)

                                new_mol = Chem.MolFromSmiles(new_smi, sanitize=False)
                                new_mol = remove_atom_map(new_mol)

                                for atom in new_mol.GetAtoms():
                                    atom.SetIsotope(max_atom_map)
                                    max_atom_map+=1

                                smiles_w_isotope = '.'.join([node['smiles_w_isotope'], Chem.MolToSmiles(new_mol)])
                                smiles_w_mapping = Chem.MolToSmiles(isotope_to_atommap(Chem.MolFromSmiles(smiles_w_isotope, sanitize=False)))
                                mol = Chem.MolFromSmiles(smiles_w_isotope, sanitize=False)

                                node['smiles_w_isotope'] = smiles_w_isotope
                                node['smiles_w_mapping'] = smiles_w_mapping
                                node['smiles'] = get_plain_smiles(Chem.MolFromSmiles(smiles_w_mapping, sanitize=False))
                                node['mol'] = [Chem.MolFromSmiles(smi, sanitize=False) for smi in smiles_w_isotope.split('.')]
                                num_match += 1
                                # print('nnew smi', Chem.MolToSmiles(new_mol))
                                possible_reactant_list.append(new_mol)
                                # print('possible_reactant_list', possible_reactant_list)
                                templ_mol_pair[pat] = possible_reactant_list
                                break
            
                if num_match == len(patterns):
                    template_reactant_dict[templ] = templ_mol_pair

        temple_combi_dict = {}
        # print('template_reactant_dict', template_reactant_dict)
        for templ, reactant_dict in template_reactant_dict.items():
            reactant_dict = self.deduplicate_reactants(reactant_dict)
            # print('reactant_dict', reactant_dict)
            # rmol_combinations = [list(combination) for combination in itertools.product(*reactant_dict.values()) if
            #                     len(combination) == len(reactant_dict)]
            rmol_combinations = [
                                list(permutation)
                                for combination in itertools.product(*reactant_dict.values())
                                for permutation in itertools.permutations(combination)
                                ]

            temple_combi_dict[templ] = rmol_combinations

        return temple_combi_dict, node
    
    def deduplicate_reactants(self, reactant_dict):
        from collections import defaultdict
        grouped_keys = defaultdict(list)
        for key, value_list in reactant_dict.items():
            grouped_keys[tuple(value_list)].append(key)
        result_dict = {}
        for value_list, keys in grouped_keys.items():
            unique_values = list(value_list)  # Unique values (n)
            m = len(keys)  # Number of keys
            n = len(unique_values)  # Number of unique values
            
            # Step 3: Assign unique values to keys
            for i, key in enumerate(keys[:min(m, n)]):
                result_dict[key] = [unique_values[i]]
            
            # Step 4: Distribute remaining values evenly
            remaining_keys = keys[min(m, n):]
            remaining_values = unique_values[min(m, n):]
            if remaining_keys:
                per_key = len(remaining_values) // len(remaining_keys)
                extra = len(remaining_values) % len(remaining_keys)
                
                idx = 0
                for key in remaining_keys:
                    end_idx = idx + per_key + (1 if extra > 0 else 0)
                    result_dict[key] = remaining_values[idx:end_idx]
                    idx = end_idx
                    if extra > 0:
                        extra -= 1

        return result_dict

    # def allow_duplicating_reactant(self, problem_temple_reactant_dict, G):
    #     """
    #     When there isn't reactant pair but some already consumed reactant,
    #     allowing duplications of reactant.
    #     This code duplicates only the molecule having the identity as reactant.
    #     """
    #     reactant_idx = set(flatten_list([value for value in problem_temple_reactant_dict.values()][0]))
    #     args = self.args

    #     new_reactant_smiles = []

    #     for rmol in reactant_idx:
    #         r_identity = G.nodes[rmol]['mol_node'].identity
    #         if r_identity == 'reactant':
    #             new_reactant_smiles.append(G.nodes[rmol]['mol_node'].smiles)

    #     return new_reactant_smiles

    def template_enumeration(self, node):
        """Enumerate templates to balance proton and allow unimolecular reaction"""
   
        if self.args.proton:
            self.proton_balanced_reaction(node)

        if self.args.uni_rxn:
            self.uni_molecular_reaction()

        if self.args.explicit_H:
            self.get_explcit_H_template()

    def proton_balanced_reaction(self, G):
        pKa_list = self.pKa_list
        templ_list = self.template_list
        if pKa_list == [None]*len(pKa_list):
            # if self.args.verbosity:
            #     logging.info('Proton balancing is not needed for this step')
            return

        new_rxn_template = []
        for pKa, templ in zip(pKa_list, templ_list):

            if not pKa:
                new_rxn_template.append(templ)
                continue

            A = pKa.get('A')
            B = pKa.get('B')

            possible_acid_base = []
            if A:
                filtered_data = [d for d in AcidBase_lookup.Acid_base if 'A' in d['role'] and d['pKa'] <= A]
                sorted_data = sorted(filtered_data, key=lambda x: x['pKa'])
                possible_acid_base = find_acid_base(G, sorted_data, 'A')
            if B:
                filtered_data = [d for d in AcidBase_lookup.Acid_base if 'B' in d['role'] and d['pKa'] >= B]
                sorted_data = sorted(filtered_data, key=lambda x: x['pKa'], reverse=True)
                possible_acid_base = find_acid_base(G, sorted_data, 'B')
            if not possible_acid_base: continue

            for acid_base in possible_acid_base:
                # for templ in rxn_templates:
                r, p = templ.split('>>')
                new_r = '.'.join([r, acid_base[0]])
                new_p = '.'.join([p, acid_base[1]])
                new_templ='>>'.join([new_r,new_p])
                new_rxn_template.append(new_templ)

        if self.args.verbosity:
            logging.info(f'{len(new_rxn_template)} new proton balanced templates were made from {len(templ_list)} templates')

        if not new_rxn_template:
            if self.args.verbosity:
                logging.info(f'No acid base error')
            # raise NoAcidBaseError

        self.template_list = new_rxn_template

    def uni_molecular_reaction(self):
        templ_list = self.template_list

        new_templates = []
        for templ in templ_list:
            rxn = AllChem.ReactionFromSmarts(templ)
            num_reactants = rxn.GetNumReactantTemplates()

            if num_reactants > 1:
                new_templates.append(templ)
                r, p = templ.split('>>')
                r = '(' + r + ')'
                p = '(' + p + ')'
                new_temp = '>>'.join([r, p])
                new_templates.append(new_temp)
            else:
                new_templates.append(templ)

        self.template_list = new_templates

    def get_explcit_H_template(self):

        new_templ_list = []
        for templ in self.template_list:
            new_templ, _ = modify_explicit_H(templ)
            new_templ_list.append(new_templ)

        self.template_list = new_templ_list

def find_acid_base(node, filtered_list, ab_condition):
    possible_acid_base=[]
    # print(filtered_list)

    for acid_base in filtered_list:
        if ab_condition=='A':
            reactant=acid_base['Acid']
            product=acid_base['Base']
        elif ab_condition=='B':
            reactant=acid_base['Base']
            product=acid_base['Acid']

        patt = Chem.MolFromSmarts(reactant)
        patt.UpdatePropertyCache(strict=False)

        mols = node['mol']
        # print('mols in acid_base', [Chem.MolToSmiles(mol) for mol in mols])
        # print(mols)

        for mol in mols:
            mol.UpdatePropertyCache(strict=False)
            try:
                if mol and mol.GetSubstructMatch(patt) and [reactant, product] not in possible_acid_base:
                    # print(reactant)
                    possible_acid_base.append([reactant, product])
            except: continue #TODO: CO makes error with a pattern of '[O:96]=[CH:97][C:99]([NH2:98])=[O:100]'
    return possible_acid_base
