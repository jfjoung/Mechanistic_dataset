
import os
import logging
import json
import pickle
from multiprocessing import Pool, Value

# Third-party libraries
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import networkx as nx
from tqdm import tqdm

# Custom modules
from scripts.setting import Args
from scripts.generate_mech_data import reagent_matching_for_single_reaction, calling_rxn_template
from utils import template_process2
from utils.template_process2 import remove_atom_map, isotope_to_atommap, get_plain_smiles, flatten_list

# Disable RDKit warnings
RDLogger.DisableLog('rdApp.*')

# NUM = 0

class Reaction_Network:
    def __init__(self, rxn_dict, idx, args):
        self.reaction_class = rxn_dict['reaction_name']
        self.reaction_smiles = rxn_dict['reaction_smiles']
        self.reaction_condition = rxn_dict['conditions'][idx]
        self.args = args
        self.rxn_network = nx.DiGraph()
        self.pruned_graph = None
        self.product_nodes = None
        self.atom_mapping()


    def atom_mapping(self):
        G = self.rxn_network
        original_rxn = self.reaction_smiles

        reactants, reagents, psmi =  original_rxn.split('>')
        if reagents:
            rsmi = '.'.join([reactants, reagents])
        else:
            rsmi = reactants

        self.ps = Chem.SmilesParserParams()
        self.ps.sanitize = False

        rmol = Chem.MolFromSmiles(rsmi, sanitize=False)
        rmol = remove_atom_map(rmol)
        pmols = [Chem.MolFromSmiles(smi, sanitize=False) for smi in psmi.split('.')]
        # pmol = Chem.MolFromSmiles(psmi, sanitize=False)
        # pmol = remove_atom_map(pmol)
        
        if len(pmols) > 1:  #TODO: If products are multiple, it needs to pick the major product.
            pmol_dict = {}
            for p in pmols:
                pmol_dict[p] = p.GetNumHeavyAtoms()
            pmol = max(pmol_dict, key=pmol_dict.get)
            pmol = remove_atom_map(pmol)
        else: 
            pmol = pmols[0]
            pmol = remove_atom_map(pmol)

        self.pmol = pmol
        self.psmi = Chem.MolToSmiles(pmol, isomericSmiles=False)
        # print('self.psmi', self.psmi)
        

        if args.explicit_H:
            self.ps.removeHs = False
            rmol.UpdatePropertyCache(strict=False)
            rmol = Chem.AddHs(rmol, explicitOnly=False)

            hpmols = [Chem.MolFromSmiles(smi, sanitize=False) for smi in psmi.split('.')]

            if len(hpmols) > 1:  #TODO: If products are multiple, it needs to pick the major product.
                pmol_dict = {}
                for p in pmols:
                    pmol_dict[p] = p.GetNumHeavyAtoms()
                hpmol = max(pmol_dict, key=pmol_dict.get)
                hpmol = remove_atom_map(hpmol)
            else: 
                hpmol = hpmols[0]
                hpmol = remove_atom_map(hpmol)
            hpmol.UpdatePropertyCache(strict=False)
            hpmol = Chem.AddHs(hpmol, explicitOnly=False)
            self.pmol = hpmol
            self.psmi = Chem.MolToSmiles(hpmol)
            # print('self.psmi', self.psmi)

        atom_map = 1
        for atom in rmol.GetAtoms():
            atom.SetIsotope(atom_map)
            atom_map+=1
        # print(atom_map)
        # print(Chem.MolToSmiles(rmol))
        # Generate the SMILES and molecular representations
        smiles_w_isotope = Chem.MolToSmiles(rmol) #)
        # print('smiles_w_isotope', smiles_w_isotope)
        smiles_w_mapping = Chem.MolToSmiles(isotope_to_atommap(rmol)) #)
        # print('smiles_w_mapping', smiles_w_mapping)
        mol = Chem.MolFromSmiles(smiles_w_isotope, sanitize=False)

        self.node_id = 0
        G.add_node(self.node_id)
        # Add these attributes to the node in the NetworkX graph
        G.nodes[self.node_id]['smiles_w_isotope'] = smiles_w_isotope
        G.nodes[self.node_id]['smiles_w_mapping'] = smiles_w_mapping
        G.nodes[self.node_id]['smiles'] = get_plain_smiles(mol)
        G.nodes[self.node_id]['mol'] = [Chem.MolFromSmiles(smi, sanitize=False) for smi in smiles_w_isotope.split('.')]
        # print('Mol', [Chem.MolToSmiles(mol) for mol in G.nodes[self.node_id]['mol']])
        G.nodes[self.node_id]['type'] = 'mol_node'
        # self.print_graph()

    def print_graph(self):
        if not self.args.do_not_pruning and self.pruned_graph:
            G = self.pruned_graph
        else:
            G = self.rxn_network

        print("Nodes of the graph:")
        for node, data in G.nodes(data=True):
            print('node', node)
            print(data)
        print("\nEdges of the graph:")
        for edge in G.edges(data=True):
            print(edge)
        

    def set_template_dict(self, templ_dict):
        self.reagents = templ_dict['Reagent'] # This is already used in reagent_matching_for_single_reaction
        self.exclude_reagents = templ_dict['Exclude_reagent'] # This is already used in reagent_matching_for_single_reaction
        self.stage = templ_dict['Stages']
        self.length = len(self.stage)

    def combine_mol(self,mol_list):
        if len(mol_list) == 1:
            return Chem.MolToSmiles(mol_list[0])
        else:
            smiles = '.'.join([Chem.MolToSmiles(mol) for mol in mol_list])
            return smiles
        
    def run_reaction(self):
        G = self.rxn_network

        for step in range(self.length):
            frontier_nodes = [node for node in G if G.out_degree(node) == 0]
            if self.args.verbosity:
                print(step, frontier_nodes)
            if len(frontier_nodes) > self.args.num_combination:
                continue
            for node in frontier_nodes:
                elementary_step = self.stage[step]
                template = template_process2.Template_process(elementary_step, self.args)

                template.template_enumeration(G.nodes[node], G.nodes[0])
                template_reactant_dict, new_node = template.find_reactants(G.nodes[node], G.nodes[0])
                for key, value in new_node.items():
                    G.nodes[node][key] = value
                if self.args.verbosity: 
                    print_plain_smiles = []
                    if template_reactant_dict:
                        print('template_reactant_dict', template_reactant_dict)
                    # for templ, mol_pairs in template_reactant_dict.items():
                    #     print(f"Template: {templ}")
                    #     for mol_list in mol_pairs:
                    #         smiles_list = [Chem.MolToSmiles(mol) for mol in mol_list]
                    #         print("  Reactants as SMILES:", smiles_list)
                # print('smiles_w_mapping',G.nodes[node]['smiles_w_mapping'])
                # print('smiles_w_isotope',G.nodes[node]['smiles_w_isotope'])
                # print(G.nodes[node]['smiles'])

                for templ, reactants in template_reactant_dict.items():
                    # print(len(reactants))

                    if len(reactants) > self.args.num_combination:
                        continue
                        # raise ValueError('Too many combinations')
                    rxn = AllChem.ReactionFromSmarts(templ)

                    for reactant in reactants:

                        combined_reactant = self.combine_mol(reactant)
                        reactant_formula = CalcMolFormula(Chem.MolFromSmiles(combined_reactant, self.ps))

                        rmols = [mol for mol in reactant]
                        # print(rmols)
                        reactant_isotopes = set([atom.GetIsotope() for mol in reactant for atom in mol.GetAtoms() ])
                        # print('reactant_isotopes', reactant_isotopes)
                        
                        try:
                            outcomes = rxn.RunReactants(rmols)
                            # print(outcomes)

                        except Exception as e:
                            # print(e)
                            continue

                        if not outcomes:
                            continue
                        # if self.args.verbosity:
                            # printing_rsmi = [G.nodes[idx]['mol_node'].__str__() for idx in reactant] #TODO: printing the reactants
                            # logging.info(f'Reactants are {printing_rsmi}')

                        for outcome in outcomes:
                            rxn_node = {'template': templ, 'description': template.description}

                            if self.args.verbosity: 
                                logging.info(f'{templ} were used')

                            # from rdkit.Chem.rdMolDescriptors import CalcMolFormula
                            # [prod.UpdatePropertyCache(strict=False) for prod in outcome]
                            # print([CalcMolFormula(prod) for prod in outcome])
                                
                            prod_smi = '.'.join([Chem.MolToSmiles(prod) for prod in outcome]) #
                            # if self.args.verbosity:
                            #     print('prod_smi', prod_smi)

                            prod_mol = Chem.MolFromSmiles(prod_smi, self.ps)
                            prod_mol.UpdatePropertyCache(strict=False)
                            product_formula = CalcMolFormula(prod_mol)
                            product_isotopes = [atom.GetIsotope() for mol in outcome for atom in mol.GetAtoms()]
                            # print('product_isotopes', product_isotopes)

                            set_product_isotopes = set(product_isotopes)
                            if len(product_isotopes) != len(set_product_isotopes):
                                continue
                            if set_product_isotopes != reactant_isotopes:
                                continue

                            if reactant_formula != product_formula:
                                continue

                            # print('product_isotopes', product_isotopes)
                            # print('prod_smi', prod_smi)

                            smi_isotope_reactant = set(G.nodes[node]['smiles_w_isotope'].split('.'))
                            reactant_smi = set([Chem.MolToSmiles(mol) for mol in reactant]) #
                            spectator_smi = smi_isotope_reactant - set([Chem.MolToSmiles(mol) for mol in reactant])
                            # print('spectator_smi', spectator_smi)
                            prod_spec_smi = '.'.join(list(spectator_smi)+[prod_smi])
                            # print(prod_spec_smi)

                            plain_smiles = get_plain_smiles(Chem.MolFromSmiles(prod_spec_smi, sanitize=False))
                            # print('plain_smiles', plain_smiles)
                            if self.args.verbosity and plain_smiles not in print_plain_smiles: 
                                print('plain_smiles:', plain_smiles)
                                print_plain_smiles.append(plain_smiles)

                            matching_nodes = [idx for idx in G.nodes 
                                            if G.nodes[idx].get('type') == 'mol_node' and G.nodes[idx].get('smiles') == plain_smiles
                                            ]
                            
                            # print('matching_nodes', matching_nodes)
                            
                            if not matching_nodes:

                                rxn_node_count = len([node for node, data in G.nodes(data=True) if data['type'] == 'rxn_node'])
                                G.add_node(f'rxn {rxn_node_count}', rxn_node=rxn_node, type='rxn_node')
                                G.add_edge(node, f'rxn {rxn_node_count}')

                                self.node_id += 1
                                G.add_node(self.node_id)
                                G.nodes[self.node_id]['smiles_w_isotope'] = prod_spec_smi
                                G.nodes[self.node_id]['smiles_w_mapping'] = Chem.MolToSmiles(isotope_to_atommap(Chem.MolFromSmiles(prod_spec_smi, sanitize=False)))
                                G.nodes[self.node_id]['smiles'] = get_plain_smiles(isotope_to_atommap(Chem.MolFromSmiles(prod_spec_smi, sanitize=False)))
                                G.nodes[self.node_id]['mol'] = [Chem.MolFromSmiles(smi, sanitize=False) for smi in prod_spec_smi.split('.')]
                                G.nodes[self.node_id]['type'] = 'mol_node'

                                # if self.args.verbosity:
                                #     print(G.nodes[self.node_id]['smiles'])


                                G.add_edge(f'rxn {rxn_node_count}', self.node_id)

                            else:
                                for node_idx in matching_nodes:
                                    # Find the parent nodes of the `matching_node`
                                    parent_nodes = [n for n in G.predecessors(node_idx) if G.nodes[n].get('type') == 'rxn_node']

                                    grandparent_nodes = []
                                    for parent in parent_nodes:
                                        # Find the grandparent nodes (predecessors of the parent)
                                        grandparent_nodes.extend([n for n in G.predecessors(parent) if G.nodes[n].get('type') == 'mol_node'])

                                    if node not in grandparent_nodes:
                                        rxn_node_count = len([node for node, data in G.nodes(data=True) if data['type'] == 'rxn_node'])
                                        G.add_node(f'rxn {rxn_node_count}', rxn_node=rxn_node, type='rxn_node')
                                        G.add_edge(node, f'rxn {rxn_node_count}')
                                        G.add_edge(f'rxn {rxn_node_count}', node_idx)

        # self.print_graph()
                                

    def find_product(self):
        # print(self.psmi)
        G = self.rxn_network
        frontier_nodes = [node for node in G if G.out_degree(node) == 0]
        product_nodes = []
        self.ps.sanitize = True

        psmi =  Chem.MolToSmiles(Chem.MolFromSmiles(self.psmi,self.ps), isomericSmiles=False,canonical=True)
        # print('psmi', psmi)
        for node_id in frontier_nodes:
            # print(node_id)
            # print(G.nodes[node_id]['smiles'])
            try:
                smi = [Chem.MolToSmiles(Chem.MolFromSmiles(smi,self.ps), isomericSmiles=False, canonical=True) for smi in G.nodes[node_id]['smiles'].split('.')]
            except Exception as e:
                # print(e)
                continue

            if psmi in smi:
                product_nodes.append(node_id)

        # patt = Chem.MolToSmiles(Chem.MolFromSmiles(self.psmi,sanitize=False))
        # patt.UpdatePropertyCache(strict=False)

        # for node_id in frontier_nodes:
        #     for mol in G.nodes[node_id]['mol']:
        #         mol.UpdatePropertyCache(strict=False)
        #         if mol.GetSubstructMatch(patt):
        #             if self.args.verbosity:
        #                 logging.info('Product is formed')
        #             product_nodes.append(node_id)

            # if self.psmi in G.nodes[node_id]['smiles']:
            #     product_nodes.append(node_id)
        # print('product_nodes', product_nodes)
        
        if not self.args.do_not_pruning:
            self.product_nodes = product_nodes
            if self.product_nodes:
                self.pruning()

    def pruning(self):
        G = self.rxn_network
        # print('self.product_nodes', self.product_nodes)

        pruned_graph = []
        for pnode in self.product_nodes:
            path_node = []
            try:
                paths = nx.all_shortest_paths(G, source=0, target=pnode)
            except nx.NetworkXNoPath:
                continue
            for path in paths:
                for nn in path:
                    path_node.append(nn)
            pruned_graph.append(nx.DiGraph(G.subgraph(path_node)))
        
        # print(pruned_graph)
        self.pruned_graphs = pruned_graph
        # print('pruned')
        # self.print_graph()
        self.check_stoichiometry()

    def check_stoichiometry(self):
        # G = 

        # self.print_graph()

        for G in self.pruned_graphs:

            reactant_mapping_numbers = set()
            reactant_smiles = G.nodes[0]['smiles_w_mapping']  #reactant node
            reactant_mapping_numbers.update(extract_atom_mapping_numbers(reactant_smiles))

            mol_node_indices = [idx for idx, data in G.nodes(data=True) if data.get('type') == 'mol_node']

            visited_node = []
            for node_idx in mol_node_indices:

                node_smi_list = G.nodes[node_idx]['smiles_w_isotope'].split('.')

                for node_smi in node_smi_list:
                    node_smi_atom_map =  Chem.MolToSmiles(isotope_to_atommap(Chem.MolFromSmiles(node_smi, sanitize=False)))
                    atom_map = extract_atom_mapping_numbers(node_smi_atom_map)
                    if atom_map - reactant_mapping_numbers:
                        for v_node in visited_node:
                            G.nodes[v_node]['smiles_w_isotope'] = '.'.join([G.nodes[v_node]['smiles_w_isotope'], node_smi])
                            G.nodes[v_node]['smiles_w_mapping'] = Chem.MolToSmiles(isotope_to_atommap(Chem.MolFromSmiles(G.nodes[v_node]['smiles_w_isotope'], sanitize=False)))
                            G.nodes[v_node]['smiles'] = get_plain_smiles(isotope_to_atommap(Chem.MolFromSmiles(G.nodes[v_node]['smiles_w_isotope'], sanitize=False)))
                            G.nodes[v_node]['mol'] = [Chem.MolFromSmiles(smi, sanitize=False) for smi in G.nodes[v_node]['smiles_w_isotope'].split('.')]
                        reactant_mapping_numbers.update(atom_map)

                visited_node.append(node_idx)
            
        # self.print_graph()

    def get_rxn(self):

        global NUM
        
        elem_rxn = []
        root_dict ={}
        # NUM += 1


        for G in self.pruned_graphs[:self.args.max_pathways]:
            with NUM.get_lock():
                if G.nodes[0]['smiles'] in root_dict:
                    num_rxn = root_dict[G.nodes[0]['smiles']]
                else:
                    NUM.value += 1
                    root_dict[G.nodes[0]['smiles']] = NUM.value
                    num_rxn = NUM.value
            rxn_node_indices = sorted(idx for idx, data in G.nodes(data=True) if data.get('type') == 'rxn_node')
            if len(rxn_node_indices) >= self.args.num_reaction_node:
                continue
            for r_node in rxn_node_indices:
                # Get parent nodes (predecessors)
                parent_node = list(G.predecessors(r_node))[0]
                # Get child nodes (successors)
                child_node = list(G.successors(r_node))[0]
                reactant_smi = G.nodes[parent_node]['smiles_w_mapping']
                product_smi = G.nodes[child_node]['smiles_w_mapping']
                rxn_smi = f'{reactant_smi}>>{product_smi}'

                if args.rxn_numbering:
                    des = G.nodes[r_node]['rxn_node']['description']
                    rxn_smi = f'{rxn_smi}|{self.reaction_class}|{self.reaction_condition}|{des}|{num_rxn}'

                elem_rxn.append(rxn_smi)

            if self.args.end:
                for pnode in self.product_nodes:
                    try:
                        reactant_smi = G.nodes[pnode]['smiles_w_mapping']
                        rxn_smi = f'{reactant_smi}>>{reactant_smi}'

                        if args.rxn_numbering:
                            des = 'End of reaction'
                            rxn_smi = f'{rxn_smi}|{self.reaction_class}|{self.reaction_condition}|{des}|{num_rxn}'

                        elem_rxn.append(rxn_smi)
                    except:
                        continue

        return elem_rxn
    

def extract_atom_mapping_numbers(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    return {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}




def generate_elementary_reaction(input):
    """
    It makes elementary reactions for a single overall reactions.
    It takes a line of 'reactants>>products NameRXN_class'
    It returns two dictionaries for the reaction networks and statistics
    """
    line, args = input

    rxn = line.split()[0]
    if len(line.split()[1:]) == 1:
        label = line.split()[-1]
    else:
        label = ' '.join(line.split()[1:])

    # When reaction class is specified, then corresponding reactions will only be considered.
    if args.rxn_class and args.rxn_class not in label:
        return None
    if args.exclude_rxn and label in args.exclude_rxn:
        return None

    rxn_dict = {'reaction_name': label,
                'reaction_smiles': rxn, }

    # try:
    rxn_dict = reagent_matching_for_single_reaction(rxn_dict, label)
    # print(rxn_dict)
    if not rxn_dict['conditions']:
        statistics = {'No templates' : 1}
        return rxn, []
    results = get_mechanistic_network(rxn_dict, args)
    # except Exception as e:
    #     print(e)
    #     if args.verbosity:
    #         logging.info(f'{e}')
    #     return rxn, []

    return rxn, results


def get_mechanistic_network(rxn, args):
    conditions = calling_rxn_template(rxn)
    if args.verbosity:
        logging.info('Generating mechanism for: ' + rxn['reaction_smiles'])
    elem_rxns = []
    for i, cond in enumerate(conditions):
        # try:
        reaction = Reaction_Network(rxn, i, args)
        reaction.set_template_dict(cond)
        reaction.run_reaction()
        reaction.find_product()
        if reaction.product_nodes:
            elem_rxn = reaction.get_rxn()
            elem_rxns.append(elem_rxn)

        # except Exception as e: 
        # #     # print(e)
        #     continue

    return flatten_list(elem_rxns)

def generate_mechdata_multiprocess(args):
    '''
    Input format is a string of 'reaction_smiles NameRXN_name'
    '''

    # Make directory for saving file
    save_dir_path=os.path.dirname(args.save)
    os.makedirs(save_dir_path, exist_ok=True)

    # def initializer():
    #     global NUM

    def init(arg):
        global NUM
        NUM = arg

    NUM = Value('i', 0)
    
    p = Pool(args.process, initializer=init, initargs=(NUM, ))
    num_generated_rxn=0

    with open(args.data, 'r') as file:
        lines = file.readlines()
        iterables = [(line, args) for line in lines]

    if args.stat:
        statistics = {}

    base_file_root, _ = os.path.splitext(args.save)
    debug_file_path = f"{base_file_root}_debug.txt"

    if args.all_info:
        writing_format = 'wb'
    else: writing_format = 'w'

    with open(args.save, writing_format) as fout, open(debug_file_path, 'w') as file_debug:
        for results in tqdm(p.imap_unordered(generate_elementary_reaction, iterables), total=len(lines)):
            if results is not None:
                rxn, elem_rxns = results
                if elem_rxns:
                    # num_used_rxn += 1
                    if args.all_info:
                        pickle.dump(elem_rxns, fout, protocol=pickle.HIGHEST_PROTOCOL)
                    else:
                        num_generated_rxn += len(elem_rxns)
                        if args.rxn_numbering:
                            # fout.write(f'|{num_used_rxn}\n'.join(elem_rxns) + f'|{num_used_rxn}\n')
                            fout.write('\n'.join(elem_rxns) + '\n')
                        else:
                            fout.write('\n'.join(elem_rxns) + '\n')

                else:
                    # file_debug.write(f"{rxn} {stat}\n")
                    continue
                    
                # if args.stat:
                #     merge_dicts(statistics, stat)

        # if args.stat:
        #     base_file_root, _ = os.path.splitext(args.save)
        #     stat_file_path = f"{base_file_root}_statistics.json"
        #     with open(stat_file_path, 'w') as file:
        #         json.dump(statistics, file, indent=4)


    # if not args.debug:
    #     os.remove(debug_file_path)

    # logging.info(f'In total, {num_used_rxn} overall reactions were utilized')
    # logging.info(f'{num_generated_rxn} elementary steps were generated')
    # logging.info(f'Mechanistic dataset generation is done')


if __name__ == '__main__':
    args=Args()
    generate_mechdata_multiprocess(args)