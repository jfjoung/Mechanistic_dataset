from utils import AcidBase_lookup
import itertools
import networkx as nx
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from utils.exceptions import *
from utils.explicit_H import  modify_explicit_H
RDLogger.DisableLog('rdApp.*')

def flatten_list(lst):
    flattened = []
    for item in lst:
        if isinstance(item, list):
            flattened.extend(flatten_list(item))
        else:
            flattened.append(item)
    return flattened


def has_duplicate_atom(atom_map_list):
    for sublist in atom_map_list:
        if len(sublist) != len(set(sublist)):
            return True
    return False

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

    def __str__(self):
        return f'Templates of elementary step for {self.description}'

    def check_num_templates(self):
        if len(self.template_list) > self.args.max_num_temp:
            raise TooManyTemplates

    def find_reactants(self, G):

        # self.template_enumeration(G)
        self.check_num_templates()

        template_reactant_dict = {} #Get every reactant for each template
        for idx, templ in enumerate(self.template_list):
            rxn = AllChem.ReactionFromSmarts(templ)
            patterns = [rmol for rmol in rxn.GetReactants()] #Get reactants as described in the templates
            # mol_in_flask = rxn_flask.flask

            templ_mol_pair = {}  #Save the pattern - molecule pair
            num_match = 0

            for pat in patterns:
                possible_reactant_list = []  #Save the matched reactant
                for node, data in G.nodes(data=True):
                    if data['type'] == 'mol_node':
                        mol_node = data['mol_node']
                        mol = mol_node.mol

                        if mol and mol.GetSubstructMatch(pat): # and mol_node.smiles not in possible_reactant_smiles_list:
                            possible_reactant_list.append(node)
                            # print('Mol :', Chem.MolToSmiles(mol), 'Pat :', Chem.MolToSmarts(pat))

                if possible_reactant_list:
                    num_match += 1
                    templ_mol_pair[pat] = possible_reactant_list

            if num_match == len(patterns):
                template_reactant_dict[templ] = templ_mol_pair

        if not template_reactant_dict:
            raise NoReactantError
        else:
            template_reactant_dict, problem_temple_reactant_dict = self.check_reaction_history(template_reactant_dict, G)
            return template_reactant_dict, problem_temple_reactant_dict

    def check_reaction_history(self, template_reactant_dict, G):
        """
        Prepare reactant pair if there are two reactants in one template,
        Check if the same atom mapping number is in both reactants,
        Check if one of the reactant is a product of another.
        """

        args = self.args
        clean_temple_reactant_dict = {}
        problem_temple_reactant_dict = {} #Save problematic ones when duplicating the reactant is allowed.

        for templ, reactant_dict in template_reactant_dict.items():
            'Make all combinations'
            rmol_combinations = [list(combination) for combination in itertools.product(*reactant_dict.values()) if
                            len(combination) == len(reactant_dict)]
            clean_combinations = []
            problem_combinations = []
            possible_reactant_smiles_list = []
            # print('rmol_combinations', rmol_combinations)

            for combi in rmol_combinations:
                has_problem = False
                'For each combination, check if it is unrealistic'
                if len(combi) > 1:
                    'Check if they are connected'
                    for combi_for_connectivity in itertools.permutations(combi, 2):
                        combi1, combi2 = combi_for_connectivity
                        # print('Combinations', G.nodes[combi1]['mol_node'], G.nodes[combi2]['mol_node'])
                        try:
                            paths_between_two = [path for path in nx.all_shortest_paths(G, source=combi1, target=combi2)]
                            if paths_between_two and combi not in problem_combinations:
                                problem_combinations.append(combi)
                                has_problem = True
                                break
                        except nx.NetworkXNoPath:
                            # print('No path')
                            continue

                'Check if they share the same mapping number'
                reactant_atommap = [list(i) for i in [G.nodes[idx]['mol_node'].atom_mapping for idx in combi]]

                if has_duplicate_atom(reactant_atommap): #Filter if there are two or more the same atom map numbers present in the reactants
                    problem_combinations.append(combi)
                    has_problem = True
                    continue

                rsmiles_set = set([G.nodes[combi1]['mol_node'].smiles for combi1 in combi])
                if not has_problem and rsmiles_set not in possible_reactant_smiles_list:
                    clean_combinations.append(combi)
                    possible_reactant_smiles_list.append(rsmiles_set)

            if clean_combinations:
                clean_temple_reactant_dict[templ] = clean_combinations

            if problem_combinations and args.stoichiometry:
                problem_temple_reactant_dict[templ] = problem_combinations
        return clean_temple_reactant_dict, problem_temple_reactant_dict

    def allow_duplicating_reactant(self, problem_temple_reactant_dict, G):
        """
        When there isn't reactant pair but some already consumed reactant,
        allowing duplications of reactant.
        This code duplicates only the molecule having the identity as reactant.
        """
        reactant_idx = set(flatten_list([value for value in problem_temple_reactant_dict.values()][0]))
        args = self.args

        new_reactant_smiles = []

        for rmol in reactant_idx:
            r_identity = G.nodes[rmol]['mol_node'].identity
            if r_identity == 'reactant':
                new_reactant_smiles.append(G.nodes[rmol]['mol_node'].smiles)

        return new_reactant_smiles

    def template_enumeration(self, G):
        """Enumerate templates to balance proton and allow unimolecular reaction"""

        if self.args.proton:
            self.proton_balanced_reaction(G)

        if self.args.uni_rxn:
            self.uni_molecular_reaction()

        if self.args.explicit_H:
            self.get_explcit_H_template()

    def proton_balanced_reaction(self, G):
        pKa_list = self.pKa_list
        templ_list = self.template_list
        if pKa_list == [None]*len(pKa_list):
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
        if not new_rxn_template:
            raise NoAcidBaseError

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

def find_acid_base(G, filtered_list, ab_condition):
    possible_acid_base=[]

    for acid_base in filtered_list:
        if ab_condition=='A':
            reactant=acid_base['Acid']
            product=acid_base['Base']
        elif ab_condition=='B':
            reactant=acid_base['Base']
            product=acid_base['Acid']

        patt = Chem.MolFromSmarts(reactant)
        patt.UpdatePropertyCache(strict=False)
        mols = [data['mol_node'].mol for node, data in G.nodes(data=True) if data['type'] == 'mol_node']
        for mol in mols:
            mol.UpdatePropertyCache(strict=False)
            try:
                if mol and mol.GetSubstructMatch(patt) and [reactant, product] not in possible_acid_base:
                    possible_acid_base.append([reactant, product])
            except: continue #TODO: CO makes error with a pattern of '[O:96]=[CH:97][C:99]([NH2:98])=[O:100]'
    return possible_acid_base
