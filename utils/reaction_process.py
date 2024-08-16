from utils import molecule_process,template_process
from utils.exceptions import *
from templates import Reaction_templates
import networkx as nx
import itertools
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
import logging
RDLogger.DisableLog('rdApp.*')

def flatten_list(lst):
    flattened = []
    for item in lst:
        if isinstance(item, list):
            flattened.extend(flatten_list(item))
        else:
            flattened.append(item)
    return flattened

def mols_from_smiles_list(all_smiles):
    """Given a list of smiles strings, this function creates rdkit
    molecules"""

    mols = []
    for smiles in all_smiles:
        if not smiles: continue
        mols.append(Chem.MolFromSmiles(smiles, sanitize=False))
    return mols

class Reaction_Network:
    def __init__(self, rxn_dict, idx, args):
        self.reaction_class = rxn_dict['reaction_name']
        self.reaction_smiles = rxn_dict['reaction_smiles']
        self.reaction_condition = rxn_dict['conditions'][idx]
        self.args = args
        self.recorded_product = False
        self.atom_map = 1
        self.rxn_network = nx.DiGraph()
        self.set_initial_network()

    def __str__(self):
        return f'The {self.reaction_class} reaction network under the {self.reaction_condition} condition'

    def set_template_dict(self, templ_dict):
        self.reagents = templ_dict['Reagent']
        self.exclude_reagents = templ_dict['Exclude_reagent']
        self.stage = templ_dict['Stages']
        self.length = len(self.stage)

    def set_initial_network(self):
        G = self.rxn_network
        rxn_smi = self.reaction_smiles

        reactants, agents, pmol = [mols_from_smiles_list(x) for x in
                                   [mols.split('.') for mols in rxn_smi.split('>')]]
        rmol = reactants + agents

        for i, _mol in enumerate(rmol):
            _mol = Chem.MolFromSmiles(Chem.MolToSmiles(_mol)) #To canonicalize the molecule
            mol_node = molecule_process.Molecule_Node(_mol, self.args)
            self.atom_map = mol_node.add_reactant(self.atom_map)
            mol_node.index = i
            G.add_node(i, mol_node=mol_node, type='mol_node')
        
        if len(pmol) > 1:
            pmol_dict = {}
            for p in pmol:
                pmol_dict[p] = p.GetNumHeavyAtoms()
            pmol = max(pmol_dict, key=pmol_dict.get)
        else: pmol = pmol[0]

        pmol = Chem.MolFromSmiles(Chem.MolToSmiles(pmol))  # TODO: If products are multiple, it needs to pick the major product.
        pmol_node = molecule_process.Molecule_Node(pmol, self.args)
        pmol_node.add_product()
        self.recorded_product = pmol_node

    def print_graph(self):
        G = self.rxn_network
        print("Nodes of the graph:")
        for node, data in G.nodes(data=True):
            if data['type'] == 'mol_node':
                print(f"{node}: {data['mol_node'].smiles_w_mapping}")
            elif data['type'] == 'rxn_node':
                print(f"{node}: {data['rxn_node']['description']}")
        print("\nEdges of the graph:")
        for edge in G.edges(data=True):
            print(edge)

    def run_reaction(self, step):
        elementary_step = self.stage[step]

        template = template_process.Template_process(elementary_step, self.args)
        template.template_enumeration(self.rxn_network)
        templ_reactant_dict, problem_dict = template.find_reactants(self.rxn_network)

        if self.args.stoichiometry and problem_dict and not templ_reactant_dict:
            new_reactant_smiles = template.allow_duplicating_reactant(problem_dict, self.rxn_network)
            self.add_duplicated_reactant(new_reactant_smiles)
            templ_reactant_dict, _ = template.find_reactants(self.rxn_network)

        if not templ_reactant_dict:
            raise NoReactantError

        G = self.rxn_network
        for templ, reactants in templ_reactant_dict.items():
            if len(reactants) > self.args.num_combination:
                raise ValueError('Too many combinations')
            rxn = AllChem.ReactionFromSmarts(templ)
            for reactant in reactants:
                rmols = [G.nodes[idx]['mol_node'].mol for idx in reactant]
                try:
                    outcomes = rxn.RunReactants(rmols)
                except Exception as e:
                    # print(e)
                    continue

                if not outcomes:
                    continue
                if self.args.verbosity:
                    printing_rsmi = [G.nodes[idx]['mol_node'].__str__() for idx in reactant]
                    logging.info(f'Reactants are {printing_rsmi}')

                for outcome in outcomes:
                    rxn_node = {'template': templ, 'description': template.description}
                    produced_molecule_idx = self.find_duplicated_outcomes(outcome)
                    # print('produced_molecule_idx', produced_molecule_idx)
                    
                    if produced_molecule_idx:
                        # Check if a connection already exists between ridx and pidx through an existing rxn_node
                        connection_exists = False
                        for rxnnode in [node for node, data in G.nodes(data=True) if data['type'] == 'rxn_node']:
                            # Check if all reactants are connected to this rxn_node
                            reactants_connected = all(rxnnode in G.successors(ridx) for ridx in reactant)
                            
                            # Check if all products are connected to this rxn_node
                            products_connected = all(pidx in G.successors(rxnnode) for pidx in produced_molecule_idx)

                            # If all reactants and products are connected through this rxn_node, set the flag
                            if reactants_connected and products_connected:
                                connection_exists = True
                                break
                        # If no existing rxn_node connects ridx and pidx, create a new rxn_node
                        if not connection_exists:
                            rxn_node_count = len([node for node, data in G.nodes(data=True) if data['type'] == 'rxn_node'])
                            
                            G.add_node(f'rxn {rxn_node_count}', rxn_node=rxn_node, type='rxn_node')
                            for ridx in reactant:
                                G.add_edge(ridx, f'rxn {rxn_node_count}')
                            for pidx in produced_molecule_idx:
                                G.add_edge(f'rxn {rxn_node_count}', pidx)

    def find_duplicated_outcomes(self, outcome):
        G = self.rxn_network
        produced_molecule_idx = []

        query_mol_nodes = [molecule_process.Molecule_Node(prod_mol, self.args) for prod_mol in outcome]
        if self.args.verbosity: 
            logging.info(f'{[mol_node.__str__() for mol_node in query_mol_nodes]} were produced')
        [mol_node.add_intermediate() for mol_node in query_mol_nodes]
        query_mol_atommap = flatten_list([mol_node.atom_mapping for mol_node in query_mol_nodes])

        if len(query_mol_atommap) != len(set(query_mol_atommap)):
            # If the outcomes have duplicated atom mapping number, discard them
            return None

        mols_in_flask = []

        for qmol in query_mol_nodes:
            mol_in_flask = [idx for idx, mol in G.nodes(data=True) if mol['type'] == 'mol_node' and mol['mol_node'].smiles == qmol.smiles]
            mols_in_flask.append(mol_in_flask)

        if len(mols_in_flask) == len(query_mol_nodes):
            combinations = list(itertools.product(*mols_in_flask))
            for combi in combinations:
                flask_atommap = flatten_list([G.nodes[idx]['mol_node'].atom_mapping for idx in combi])
                if set(query_mol_atommap) == set(flask_atommap):
                    return list(combi)  # False #No need to add

        #If there is not a set matching with a set in the flask, then examine individual chemicals
        #If chemical is in the flask, return its index
        #If not, return a new index for new chemical

        for qmol in query_mol_nodes:
            matched_mol = [idx for idx, mol in G.nodes(data=True) if
                            mol['type'] == 'mol_node' and mol['mol_node'].smiles_w_mapping == qmol.smiles_w_mapping] #Take the first one and it should be only one
            if matched_mol:
                produced_molecule_idx.append(matched_mol[0])
            else:
                pmol_idx = self.add_new_chemical(qmol)
                produced_molecule_idx.append(pmol_idx)

        return produced_molecule_idx


    def add_new_chemical(self, mol_node):
        G = self.rxn_network
        i = len([(idx, mol['mol_node']) for idx, mol in G.nodes(data=True) if mol['type'] == 'mol_node'])
        G.add_node(i, mol_node=mol_node, type='mol_node')
        mol_node.index = i
        return i

    def add_duplicated_reactant(self, smiles_List):
        for rsmi in smiles_List:
            rmol = Chem.MolFromSmiles(rsmi)
            mol_node = molecule_process.Molecule_Node(rmol, self.args)
            self.atom_map = mol_node.add_reactant(self.atom_map)
            i = len([(idx, mol['mol_node']) for idx, mol in self.rxn_network.nodes(data=True) if mol['type'] == 'mol_node'])
            self.rxn_network.add_node(i, mol_node=mol_node, type='mol_node')

    def find_product(self):
        '''
        When the reaction is over, check if the product is present in the flask
        If there is a product in the flask, change its identity from intermediate to product
        '''

        psmi = self.recorded_product.smiles
        G = self.rxn_network
        product_idx = False

        for idx, node in G.nodes(data=True):
            if node['type'] == 'mol_node':
                _smi = node['mol_node'].smiles
                if _smi == psmi:
                    product_idx = idx
                    break
        if product_idx:
            pmol_node = G.nodes[product_idx]['mol_node']
            pmol_node.identity = 'product'
            return True
        else:
            return False

    def is_product_formed(self):
        if self.find_product():
            self.assign_identity()
            if self.args.verbosity:
                logging.info(f'Product is produced')
            return True
        else: 
            if self.args.verbosity:
                logging.info(f'Product is NOT produced')
            return False

    def assign_identity(self):
        """If product is formed, assign the chemical identity, e.g. spectator, byproduct"""

        #If the chemical nodes are disconnected, they are spectators
        isolates = [node for node in nx.isolates(self.rxn_network)]
        for node in isolates:
            if self.rxn_network.nodes[node]['type'] == 'mol_node':
                self.rxn_network.nodes[node]['mol_node'].identity = 'spectator'
            else:
                raise ValueError('There is disconnected reaction nodes!')

        leaf_nodes = [node for node, out_degree in self.rxn_network.out_degree() if out_degree == 0]
        for node in leaf_nodes:
            if node not in isolates:
                if self.rxn_network.nodes[node]['type'] == 'mol_node':
                    if  self.rxn_network.nodes[node]['mol_node'].identity != 'product':
                        self.rxn_network.nodes[node]['mol_node'].identity = 'byproduct'
                else:
                    raise ValueError('There is disconnected reaction nodes!')


if __name__ == '__main__':
    import argparse
    def parse_arguments():
        parser = argparse.ArgumentParser("Set the arguments for mechanistic dataset generation")
        parser.add_argument("--explicit_H", help="SMILES with Explicit H",
                            type=bool, default=True)
        parser.add_argument("--rxn_class",
                            help="Specify the reaction class to work on. If not specified, it will work on all classes.",
                            type=str, default='')
        parser.add_argument("--stoichiometry", help="Duplicate the reactants when needed",
                            type=bool, default=True)
        parser.add_argument("--proton", help="Allow the proton-balanced reactions",
                            type=bool, default=True)
        parser.add_argument("--uni_rxn", help="Allow the unimolecular reactions",
                            type=bool, default=True)

        parser.add_argument("--max_num_temp",
                            help="The maximum number of templates allowed to consider. If uni_rxn or proton is true, the combinatorial explosion could occur",
                            type=int, default=50)
        return parser.parse_args()


    def get_mechanistic_network(rxn, conditions,  args):
        for i, cond in enumerate(conditions):
            reaction = Reaction_Network(rxn, i, args)
            reaction.set_template_dict(cond)
            for step in range(reaction.length):
                try:
                    reaction.run_reaction(step)
                except NoReactantError as e:
                    print('No reactant')
                    continue
                except NoAcidBaseError as e:
                    print('No acid base')
                    continue
            reaction.print_graph()
            print(reaction.is_product_formed())


    args = parse_arguments()
    #Cl[Pd]Cl
    line = 'BrC1=C(Br)C=CC=C1C.OB(C)O.[Cl-]>Cl[Pd]Cl.N.P(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4.P(C5=CC=CC=C5)(C6=CC=CC=C6)C7=CC=CC=C7.O>BrC8=C(C)C=CC=C8C Bromo Suzuki coupling'
    rxn = line.split()[0]
    if len(line.split()[1:]) == 1:
        label = line.split()[-1]
    else:
        label = ' '.join(line.split()[1:])

    rxn_dict = {'reaction_name': label, 'reaction_smiles': rxn, 'conditions': ['Reaction with PdX2, PPh3, and water']}
    conditions = [{'Reagent': ['[Cl,Br,I,O][Pd;+0][Cl,Br,I,O]', '[P]', '[O;H1,H2]'], 'Exclude_reagent': None, 'Stages': {0: {'Templates': ['[Cl,Br,I,O:2]-[Pd;+0:1](-[Cl,Br,I,O:3])(-[P;D4;+0:4])(-[P;D4;+0:5])>>[Cl,Br,I,O:2]-[Pd;+0:1](-[Cl,Br,I,O:3])-[P;D3;+0:4].[P;D3;+0:5]'], 'pKa': [None], 'Description': 'Palladium-ligand dissociation'}, 1: {'Templates': ['[Cl,Br,I,O:2]-[Pd;+0:1]-[Cl,Br,I,O:3].[P;D3:4]>>[Cl,Br,I,O:2]-[Pd;+0:1]-[P;H0;+1:4].[Cl,Br,I,O;-1:3]'], 'pKa': [None], 'Description': 'Palladium-ligand formation'}, 2: {'Templates': ['[OH2:1]>>[OH1;-1:1]'], 'pKa': [{'B': 14}], 'Description': 'Hydroxide ion preparation'}, 3: {'Templates': ['[Cl,Br,I,O:2]-[Pd;+0:1]-[P;H0;+1:4].[OH1;-1:5]>>[Cl,Br,I,O:2]-[Pd;+0:1]-[P;H0;+0:4]-[OH1;+0:5]', '[Cl,Br,I,O:2]-[Pd;+0:1]-[P;H0;+1:4].[OH2;+0:5]>>[Cl,Br,I,O:2]-[Pd;+0:1]-[P;H0;+0:4]-[OH2;+1:5]'], 'pKa': [None, None], 'Description': 'Hydroxide-Phosphine reaction'}, 4: {'Templates': ['[Cl,Br,I,O:2]-[Pd;+0:1]-[P;H0;+0:4]-[OH2;+1:5]>>[Cl,Br,I,O:2]-[Pd;+0:1]-[P;H0;+0:4]-[OH1;+0:5]'], 'pKa': [{'B': -1.7}], 'Description': 'Deprotonation'}, 5: {'Templates': ['[Cl,Br,I,O:2]-[Pd;H0;+0:1]-[P;H0;+0:4]-[OH1;+0:5]>>[Cl,Br,I,O:2]-[Pd;H1;+0:1].[P;H0;+0:4]=[OH0;+0:5]'], 'pKa': [None], 'Description': 'Palladium reduction'}, 6: {'Templates': ['[Cl,Br,I,O:2]-[Pd;H1;+0:1]>>[Cl,Br,I,O;H1:2].[Pd;H0;+0:1]'], 'pKa': [None], 'Description': 'Ligand dissociation - reductive elimination type'}, 7: {'Templates': ['[Cl,Br,I,O$([O]-S):7]-[#6;+0:5].[Pd&!$([Pd]-I)&!$([Pd]-Br)&!$([Pd]-Cl):6]>>[Cl,Br,I,O$([O]-S):7]-[Pd:6]-[#6;+0:5]'], 'pKa': [None], 'Description': 'Oxidative addition'}, 8: {'Templates': ['[#8;+0;H2:8].[Cl,Br,I,O$([O]-S):7]-[Pd:6]-[#6;+0:5]>>[Cl,Br,I,O$([O]-S);H1;+0:7].[#8;H1;+0:8]-[Pd:6]-[#6;+0:5]', '[#6&!$([#6]=O):10][#8;H1:8].[Cl,Br,I,O$([O]-S):7]-[Pd:6]-[#6;+0:5]>>[Cl,Br,I,O$([O]-S);H1;+0:7].[#6:10][#8:8]-[Pd:6]-[#6;+0:5]', '[Li,Na,K,Cs:10][#8:8].[Cl,Br,I,O$([O]-S):7]-[Pd:6]-[#6;+0:5]>>[Cl,Br,I,O$([O]-S);-1:7].[Li,Na,K,Cs;+1:10].[#8:8]-[Pd:6]-[#6;+0:5]', '[#6&!$([#6]=O):9]-[#8;-1:8].[Cl,Br,I,O$([O]-S):7]-[Pd:6]-[#6;+0:5]>>[Cl,Br,I,O$([O]-S);-1:7].[#6:9]-[#8;+0:8]-[Pd:6]-[#6;+0:5]'], 'pKa': [None, None, None, None], 'Description': 'Halide leaving'}, 9: {'Templates': ['[#5:2]-[#6;+0:4].[#8:8]-[Pd:6]-[#6;+0:5]>>[#6;+0:4]-[Pd:6]-[#6;+0:5].[#5:2]-[#8;+0:8]', '[#7;H3;+0:1].[#8;H1:8]-[Pd:6]-[#6;+0:5]>>[#7;H2;+0:1]-[Pd:6]-[#6;+0:5].[#8;H2;+0:8]', '[#7;H2;+0:1].[#8;H1:8]-[Pd:6]-[#6;+0:5]>>[#7;H1;+0:1]-[Pd:6]-[#6;+0:5].[#8;H2;+0:8]', '[#7;H1;+0:1].[#8;H1:8]-[Pd:6]-[#6;+0:5]>>[#7;H0;+0:1]-[Pd:6]-[#6;+0:5].[#8;H2;+0:8]', '[#7;H3;+0:1].[#8;H0:8]-[Pd:6]-[#6;+0:5]>>[#7;H2;+0:1]-[Pd:6]-[#6;+0:5].[#8;H1;+0:8]', '[#7;H2;+0:1].[#8;H0:8]-[Pd:6]-[#6;+0:5]>>[#7;H1;+0:1]-[Pd:6]-[#6;+0:5].[#8;H1;+0:8]', '[#7;H1;+0:1].[#8;H0:8]-[Pd:6]-[#6;+0:5]>>[#7;H0;+0:1]-[Pd:6]-[#6;+0:5].[#8;H1;+0:8]', '[*:2]-[Pd:1]-[Cl,Br,I:3].[B;+0:4]-[B;+0:5]>>[*:2]-[Pd:1]-[B;+0:4].[Cl,Br,I;+0:3]-[B;+0:5]', '[*:2]-[Pd:1]-[O;D2:3]-[S:6].[B;+0:4]-[B;+0:5]>>[*:2]-[Pd:1]-[B;+0:4].[B;+0:5]-[O;D2:3]-[S:6]'], 'pKa': [None, None, None, None, None, None, None, None, None], 'Description': 'Transmetalation'}, 10: {'Templates': ['[#5,#6,#7;+0:4]-[Pd:6]-[#6;+0:5]>>[#5,#6,#7;+0:4]-[#6;+0:5].[Pd:6]'], 'pKa': [None], 'Description': 'Reductive elimination'}}}]

    get_mechanistic_network(rxn_dict,conditions, args)
