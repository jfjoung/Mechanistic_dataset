from utils import reaction_process
from utils.exceptions import *
from utils.reaction_process import flatten_list
from templates import Reaction_templates
from collections import defaultdict
from itertools import product
import networkx as nx
from rdkit import Chem, RDLogger
import copy
RDLogger.DisableLog('rdApp.*')


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
                # common_nodes_list_sorted = tuple(sorted(common_nodes_list))
                ints = [x for x in common_nodes_list if isinstance(x, int)]
                strs = [x for x in common_nodes_list if isinstance(x, str)]
                ints_sorted = sorted(ints)
                strs_sorted = sorted(strs)
                common_nodes_list_sorted = tuple(ints_sorted + strs_sorted)
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

class Get_Reactions:
    def __init__(self, reaction):
        self.rxn_network = reaction.rxn_network
        self.args = reaction.args
        self.reaction_class = reaction.reaction_class
        self.reaction_smiles = reaction.reaction_smiles
        self.reaction_condition = reaction.reaction_condition
        self.reaction_node = []
        self.reactant_node = []
        self.product_node = []
        self.byproduct_node = []
        self.intermediate_node = []
        self.spectator_node = []
        self.pruned_graph = None
        self.reaction_info = None
        self.rxn_smi = []
        self.find_chemical_nodes()

    def find_chemical_nodes(self):
        if self.pruned_graph:
            G = self.pruned_graph
            self.reaction_node = []
            self.reactant_node = []
            self.product_node = []
            self.byproduct_node = []
            self.intermediate_node = []
            self.spectator_node = []

        else:
            G = self.rxn_network

        for idx, node in G.nodes.data():
            if node['type'] == 'rxn_node':
                self.reaction_node.append(idx)
            else:
                identity = node['mol_node'].identity
                if identity == 'reactant':
                    self.reactant_node.append(idx)
                elif identity == 'product':
                    self.product_node.append(idx)
                elif identity == 'byproduct':
                    self.byproduct_node.append(idx)
                elif identity == 'intermediate':
                    self.intermediate_node.append(idx)
                elif identity == 'spectator':
                    self.spectator_node.append(idx)

    def print_graph(self):
        if not self.args.do_not_pruning and self.pruned_graph:
            G = self.pruned_graph
        else:
            G = self.rxn_network
        print("Nodes of the graph:")
        for node, data in G.nodes(data=True):
            if data['type'] == 'mol_node':
                print(f"{node}: {data['mol_node'].smiles_w_mapping}, {data['mol_node'].identity}")
            elif data['type'] == 'rxn_node':
                print(f"{node}: {data['rxn_node']['description']}")
            else:
                print(node, data)
        print("\nEdges of the graph:")
        for edge in G.edges(data=True):
            print(edge)

    def get_node(self, node_id):
        if not self.args.do_not_pruning and self.pruned_graph:
            G = self.pruned_graph
        else:
            G = self.rxn_network

        if G.nodes[node_id]['type'] == 'rxn_node':
            return G.nodes[node_id]['rxn_node']
        
        elif G.nodes[node_id]['type'] == 'mol_node':
            return G.nodes[node_id]['mol_node']

    def graph_pruning(self):
        G = self.rxn_network
        reaction_path = []

        cycle_pathways = [i for i in nx.simple_cycles(G)]
        shared_cycles = find_shared_nodes(cycle_pathways)
        # print(shared_cycles)

        route_to_impurity = []
        if cycle_pathways:
            """If there are cycles, 
            check if the cycle is not related to product formation, and record them in route_to_impurity
            """

            route_to_product = []
            route_length = {}
            for shared_node, cycle_list in shared_cycles.items():
                if shared_node == 'Single Cycle': pass #If there is only one cycle, ignore.
                else:
                    for cycle in cycle_list:
                        for node in cycle:
                            for pnode in self.product_node:
                                shortest_lenght = 100
                                try:
                                    paths = nx.all_shortest_paths(G, source=node, target=pnode)
                                except nx.NetworkXNoPath:
                                    continue
                                for path in paths:
                                    shortest_lenght = min(shortest_lenght, len(path))
                                route_length[node] = shortest_lenght
                node_close_to_product = min(route_length, key=route_length.get)

                for cycle in cycle_list:
                    if node_close_to_product in cycle:
                        route_to_product.append(cycle)
                        pass
                    else:
                        route_to_impurity.append(cycle)

                route_to_impurity = flatten_list(route_to_impurity)
                route_to_impurity = list(set(route_to_impurity))
                for node in shared_node:
                    if node in route_to_impurity: route_to_impurity.remove(node)

            route_to_impurity = flatten_list(route_to_impurity)
            route_to_product = flatten_list(route_to_product)

            route_to_impurity = [node for node in route_to_impurity if G.nodes[node]['type'] == 'rxn_node']
            route_to_product = [node for node in route_to_product if G.nodes[node]['type'] == 'rxn_node']

            route_to_impurity = list(set(route_to_impurity) - set(route_to_product))
            # print('route_to_impurity', route_to_impurity)
            # print('route_to_product', route_to_product)


        # Get every reaction node connecting the reactants and products
        for rnode in self.reactant_node:
            for pnode in self.product_node:
                try:
                    paths = nx.all_shortest_paths(G, source=rnode, target=pnode)
                    for path in paths:
                        rxn_path = [node for node in path if G.nodes[node]['type'] == 'rxn_node']
                        if any(node in route_to_impurity for node in rxn_path):
                            continue
                        # print(rxn_path)
                        reaction_path.append(rxn_path)
                except nx.NetworkXNoPath:
                    # print('no path')
                    continue

        reaction_path = flatten_list(reaction_path)
        reaction_path = list(set(reaction_path))
        # print('reaction_path', reaction_path)

        # Get all neighbors of the reaction nodes
        path_node = []
        for reaction_node in reaction_path:
            path_node.append(reaction_node)
            neighbor_node = [n for n in nx.all_neighbors(G, reaction_node)]
            for nn in neighbor_node:
                if nn not in path_node:
                    path_node.append(nn)

        pruned_graph = nx.DiGraph(G.subgraph(path_node))

        # self.print_graph()

        #Check if there is disconnected node in the pruned graph
        isolates = [node for node in nx.isolates(pruned_graph)]
        disconnected_intermediates = [node for node, out_degree in pruned_graph.out_degree() if out_degree == 0 and pruned_graph.nodes[node]['mol_node'].identity == 'intermediate']

        # if disconnected_intermediates or isolates:
        #     for node in disconnected_intermediates + isolates:
        #         print(G.nodes[node]['mol_node'])

        if disconnected_intermediates or isolates:
            for rxnode in self.reaction_node:
                precursor_nodes = [node for node in G.predecessors(rxnode)]

                # The case where two disconnected nodes are the reactants for the subsequent reactions
                if set(precursor_nodes).issubset(set(disconnected_intermediates + isolates)):
                    successor_nodes = [node for node in G.successors(rxnode)]
                    possible_byproduct = [node for node in G.successors(rxnode) if G.nodes[node]['mol_node'].identity in ['byproduct', 'product']]
                    if successor_nodes == possible_byproduct:
                        path_node.append(rxnode)
                        for nn in successor_nodes+precursor_nodes:
                            if nn not in path_node:
                                path_node.append(nn)
                
                # The case where one disconnected node is producing the byproduct in one step
                elif set(disconnected_intermediates + isolates) & set(precursor_nodes):
                    successor_nodes = [node for node in G.successors(rxnode)]
                    possible_byproduct = [node for node in G.successors(rxnode) if G.nodes[node]['mol_node'].identity in ['byproduct', 'product']]
                    if successor_nodes == possible_byproduct:
                        path_node.append(rxnode)
                        for nn in successor_nodes+precursor_nodes:
                            if nn not in path_node:
                                path_node.append(nn)

            pruned_graph = nx.DiGraph(G.subgraph(path_node))
            isolates = [node for node in nx.isolates(pruned_graph)]
            disconnected_intermediates = [node for node, out_degree in pruned_graph.out_degree() if
                                          out_degree == 0 and pruned_graph.nodes[node][
                                              'mol_node'].identity == 'intermediate']
            if disconnected_intermediates or isolates:
                raise ValueError(f'There is disconnected intermediate!{disconnected_intermediates, isolates}') #TODO: If there are two or more disconnected intermediates, find a route connecting them.

        

        # Those nodes were reactants, but now they are spectators
        missing_reactant_nodes = [node_id for node_id in self.reactant_node if node_id not in pruned_graph.nodes]
        for node in missing_reactant_nodes:
            node_info = G.nodes[node]
            pruned_graph.add_node(node, mol_node=node_info['mol_node'], type='mol_node')
            pruned_graph.nodes[node]['mol_node'].identity = 'spectator'

        for node in self.spectator_node:
            node_info = G.nodes[node]
            pruned_graph.add_node(node, mol_node=node_info['mol_node'], type='mol_node')

        self.pruned_graph = pruned_graph
        # self.print_graph()
        self.find_chemical_nodes()


    def get_elementary_reactions_info(self):
        """Code for obtaining the reaction information"""
        if not self.args.do_not_pruning:
            G = self.pruned_graph
        else:
            G = self.rxn_network

        args = self.args
        if len(self.reaction_node) > args.num_reaction_node:
            raise ValueError('Too many reaction nodes')

        cycle_pathways = [i for i in nx.simple_cycles(G)]

        if cycle_pathways:
            if len(cycle_pathways) > args.num_cycles:
                raise ValueError("Too many cycles.")

            """If there are 2 or more cycles, each cycle might have different reactants and byproducts. 
            As they do not affect other reactions, find all possible routes"""
            if len(cycle_pathways) > 1:
                filtered_cycles = []
                for i, cycle in enumerate(cycle_pathways):
                    filtered_cycle = [node for node in cycle if G.nodes[node]['type'] == 'rxn_node']
                    filtered_cycles.append(filtered_cycle)
            elif len(cycle_pathways) == 1:
                filtered_cycles = [node for node in cycle_pathways[0] if G.nodes[node]['type'] == 'rxn_node']
            cycle_dict=find_shared_nodes(filtered_cycles)
            reaction_nodes_outside=list(set(self.reaction_node)-set(flatten_list(filtered_cycles)))

            def create_combinations(cycle_dict):
                # Extract sequences from each group
                sequences = [seq for sequences in cycle_dict.values() for seq in sequences]
                # Create all possible combinations of one sequence from each group
                combinations = list(product(*cycle_dict.values()))
                return combinations

            combination_of_cycles=create_combinations(cycle_dict) # In case of more than 2 cycles
            reaction_paths=[flatten_list(list(combi)+reaction_nodes_outside) for combi in combination_of_cycles]
        else: reaction_paths = [self.reaction_node]

        #Sorting the reactions

        # reaction_paths = [sorted(reactions, key=lambda x: int(x.split()[1])) for reactions in reaction_paths]
        reaction_paths = [
            sorted(reactions, key=lambda x: int(x.split()[1]) if len(x.split()) > 1 else float('inf')) for reactions in
            reaction_paths]

        # reaction_paths = [sorted(reaction_paths, key=lambda x: int(x.split()[1]))]

        reaction_node_for_product = [node for node in G.predecessors(self.product_node[0])][0]


        reaction_info = {}
        for idx, rxn_path in enumerate(reaction_paths):
            reaction_info_path = {}
            for rxn_node in rxn_path:
                #Get reactants and products for this elementary step
                precursor_nodes = [node for node in G.predecessors(rxn_node)]
                # print(precursor_nodes)
                precursor_atommap = flatten_list([G.nodes[node]['mol_node'].atom_mapping for node in precursor_nodes if G.nodes[node]['type'] == 'mol_node'])
                # print(precursor_atommap)
                successor_nodes = [node for node in G.successors(rxn_node)]

                #Check any reactants are consumed before
                conmused_reactant_node = []
                for r_node in self.reactant_node:
                    for precursor in precursor_nodes:
                        if precursor in self.reactant_node:
                            conmused_reactant_node.append(precursor)
                        else:
                            is_used = False
                            if r_node != precursor:
                                try:
                                    paths_between_two = [path for path in
                                                         nx.all_shortest_paths(G, source=r_node, target=precursor)]
                                except:
                                    continue
                                for path in paths_between_two:
                                    path = [node for node in path if G.nodes[node]['type'] == 'rxn_node']
                                    not_allowed_nodes = set(path) - set(rxn_path)
                                    if not not_allowed_nodes:
                                        conmused_reactant_node.append(r_node)
                                        is_used = True
                                        break
                                if is_used: break

                not_used_reactants = set(self.reactant_node) - set(conmused_reactant_node)
                
                #Check the byproducts are formed before
                produced_byproduct_nodes = []
                for by_node in self.byproduct_node + self.product_node:
                    reaction_node_for_byproduct = [node for node in G.predecessors(by_node)][0]
                    by_node_atommap = G.nodes[by_node]['mol_node'].atom_mapping
                    # print('by_node_mapping',by_node_atommap)
                    if set(by_node_atommap) & set(precursor_atommap):
                        continue

                    if reaction_node_for_byproduct not in rxn_path:
                        #This byproduct cannot be formed in this pathway
                        continue

                    for precursor in precursor_nodes:
                        is_produced = False
                        if precursor in self.reactant_node:
                            continue
                        try:
                            paths_between_two=[path for path in nx.all_shortest_paths(G, source=reaction_node_for_byproduct, target=precursor)]
                        except Exception as e:
                            continue

                        for path in paths_between_two:
                            if reaction_node_for_product not in path:
                                is_there_reactant = any(element in self.reactant_node for element in path)
                                if is_there_reactant: continue
                                if reaction_node_for_product in path:
                                    further_reacting_intermidiate = [node for node in
                                                                     G.successors(reaction_node_for_product) if
                                                                     G.nodes[node]['mol_node'].identity == 'intermediate']
                                    if not further_reacting_intermidiate:
                                        continue
                                path=[node for node in path if G.nodes[node]['type'] == 'rxn_node']
                                not_allowed_nodes = set(path) - set(rxn_path)
                                if not not_allowed_nodes and rxn_node not in path:
                                    produced_byproduct_nodes.append(by_node)
                                    is_produced=True
                                    break
                        if is_produced: break
                # print(rxn_node, produced_byproduct_nodes, not_used_reactants)
                reaction_info_path[rxn_node] = {'reactants': precursor_nodes,
                                                'products': successor_nodes,
                                                'not used reactants': not_used_reactants,
                                                'byproducts': produced_byproduct_nodes,
                                                'spectators': self.spectator_node}
                if args.end and rxn_node == reaction_node_for_product:
                    reaction_info_path['end rxn'] = {'reactants': successor_nodes, 'products': successor_nodes,
                                                    'not used reactants': not_used_reactants,
                                                    'byproducts': produced_byproduct_nodes,
                                                     'spectators': self.spectator_node}

            reaction_info[idx] = reaction_info_path
        self.reaction_info = reaction_info
        # print(reaction_info)
        # return self.convert_to_smiles()

    def convert_to_smiles(self):
        """
        from the reaction information, convert the node index into the SMILES

        """
        args = self.args
        reaction_info = copy.deepcopy(self.reaction_info)
        if self.pruned_graph:
            G = self.pruned_graph
        else:
            G = self.rxn_network

        rxn_smiles = []
        reaction_node_for_product = [node for node in G.predecessors(self.product_node[0])][0]

        for num_path, rxn_info in reaction_info.items():
            if args.full:
                overall_reactant = []
                overall_product = []
                overall_reagent = []

            for rxn_node, mol_dict in rxn_info.items():
                reactant_side = mol_dict['reactants']
                product_side = mol_dict['products']
                reagent_side = []

                # if args.end and rxn_node == reaction_node_for_product:
                #     reactant_side = mol_dict['products']
                #     product_side = mol_dict['products']

                if args.byproduct:
                    reagent_side.append(mol_dict['byproducts'])
                    reagent_side = flatten_list(reagent_side)

                if args.spectator:
                    reagent_side.append(self.spectator_node)
                    reagent_side.append(list(mol_dict['not used reactants']))
                    reagent_side = flatten_list(reagent_side)

                if args.full:
                    overall_reactant.append([node for node in mol_dict['reactants'] if node in self.reactant_node])
                    overall_reactant = flatten_list(overall_reactant)
                    if args.spectator:
                        overall_reagent.append(self.spectator_node)
                        overall_reagent.append(list(mol_dict['not used reactants']))
                        overall_reagent = list(set(flatten_list(overall_reagent)))
                        overall_reagent = [reagent for reagent in overall_reagent if reagent not in overall_reactant]

                    if rxn_node == reaction_node_for_product:
                        overall_product.append(mol_dict['products'])
                        if args.byproduct:
                            overall_product.append(
                            [node for node in mol_dict['byproducts'] if node in self.byproduct_node])
                        overall_product = list(set(flatten_list(overall_product)))
                else:
                    if not args.reagent:
                        reactant_side += reagent_side
                        product_side += reagent_side
                        reagent_side = []
                        reactant_side = list(set(reactant_side))
                        product_side = list(set(product_side))

                    if args.plain:
                        rsmi = '.'.join([G.nodes[node]['mol_node'].smiles for node in reactant_side])
                        re_smi = '.'.join([G.nodes[node]['mol_node'].smiles for node in reagent_side])
                        psmi = '.'.join([G.nodes[node]['mol_node'].smiles for node in product_side])
                    else:
                        rsmi = '.'.join([G.nodes[node]['mol_node'].smiles_w_mapping for node in reactant_side])
                        re_smi = '.'.join([G.nodes[node]['mol_node'].smiles_w_mapping for node in reagent_side])
                        psmi = '.'.join([G.nodes[node]['mol_node'].smiles_w_mapping for node in product_side])

                    rxn_smi = f'{rsmi}>{re_smi}>{psmi}'
                    # if not args.spectator or not args.byproduct:
                    if not args.plain:
                        rxn_smi = self.remapping(rxn_smi)
                    rxn_smiles.append(rxn_smi)


        if args.full:
            if not args.reagent:
                overall_reactant += overall_reagent
                overall_product += overall_reagent
                overall_reagent = []
            if args.plain:
                rsmi = '.'.join([G.nodes[node]['mol_node'].smiles for node in overall_reactant])
                re_smi = '.'.join([G.nodes[node]['mol_node'].smiles for node in overall_reagent])
                psmi = '.'.join([G.nodes[node]['mol_node'].smiles for node in overall_product])
            else:
                rsmi = '.'.join([G.nodes[node]['mol_node'].smiles_w_mapping for node in overall_reactant])
                re_smi = '.'.join([G.nodes[node]['mol_node'].smiles_w_mapping for node in overall_reagent])
                psmi = '.'.join([G.nodes[node]['mol_node'].smiles_w_mapping for node in overall_product])
            rxn_smi = f'{rsmi}>{re_smi}>{psmi}'
            
            rxn_smiles.append(rxn_smi)

        self.rxn_smi = rxn_smiles
        
        return rxn_smiles

    def remapping(self, rxn_smi):
        if self.args.explicit_H:
            ps = Chem.SmilesParserParams()
            ps.removeHs = False
            ps.sanitize = False
        else:
            ps = Chem.SmilesParserParams()
            ps.sanitize = False

        rsmi, resmi, psmi = rxn_smi.split('>')
        rmol = Chem.MolFromSmiles(rsmi, ps)
        remol = Chem.MolFromSmiles(resmi, ps)
        pmol = Chem.MolFromSmiles(psmi, ps)

        mapping_dict = {} #old map: new map

        for idx, atom in enumerate(rmol.GetAtoms()):
            new_map = idx+1
            mapping_dict[atom.GetAtomMapNum()] = new_map
            atom.SetAtomMapNum(new_map)

        rsmi = Chem.MolToSmiles(rmol)

        for idx, atom in enumerate(remol.GetAtoms()):
            new_re_map = new_map + idx+1
            mapping_dict[atom.GetAtomMapNum()] = new_re_map
            atom.SetAtomMapNum(new_re_map)

        for atom in pmol.GetAtoms():
            atom.SetAtomMapNum(mapping_dict[atom.GetAtomMapNum()])

        resmi = Chem.MolToSmiles(remol)
        psmi = Chem.MolToSmiles(pmol)

        return f'{rsmi}>{resmi}>{psmi}'


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

if __name__ == '__main__':
    import argparse

    def str2bool(v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
    def parse_arguments():
        parser = argparse.ArgumentParser("Set the arguments for mechanistic dataset generation")

        # Arguments for applying reaction templates to get every elementary reactions
        parser.add_argument("--num_combination",
                            help="The number of possible reactions to be considered, it prevents combinatorial explosions",
                            type=int, default=30)
        parser.add_argument("--uni_rxn", help="Allow the unimolecular reactions",
                            type=str2bool, default=True)
        parser.add_argument("--proton", help="Allow the proton-balanced reactions",
                            type=str2bool, default=True)
        parser.add_argument("--max_num_temp",
                            help="The maximum number of templates allowed to consider. If uni_rxn or proton is true, the combinatorial explosion could occur",
                            type=int, default=50)
        parser.add_argument("--stoichiometry", help="Duplicate the reactants when needed",
                            type=str2bool, default=True)

        # Arguments for handling the reaction network
        parser.add_argument("--simple",
                            help="Use all simple path finding instead of shortest path during reaction network generation",
                            type=str2bool, default=False)
        parser.add_argument("--do_not_pruning",
                            help="Get full reaction network instead of the pruned network containing only paths toward product",
                            type=str2bool, default=False)

        # Arguments for preventing the combinatorial explosions in the reaction network
        parser.add_argument("--num_cycles", help="The maximum number of cycles allowed in the reaction graph",
                            type=int, default=9)
        parser.add_argument("--num_reaction_node",
                            help="The maximum number of reaction nodes allowed in the reaction graph",
                            type=int, default=50)

        # Arguments for extracting reaction SMARTS from the reaction network
        '''
        spectator: It will add spectators for all reactions. Reactants that have not yet participated in the reaction are also included. 
                   e.g. Reactants.Spectators>>Products.Spectators
        byproduct: It will add produced byproducts as spectators for downstream reaction. 
                   e.g. upstream reaction; reactants>>intermediates.byproducts
                        downstream reaction; reactants.intermediates.byproducts>>new_intermediates.byproducts
        full: It returns the overall reactions.
        end: It returns the termination reactions. e.g. products>>products
        reagent: It changes the SMARTS format from reactants.spectators>>products.spectators to reactants>spectators>products
        '''

        parser.add_argument("--byproduct", help="Add produced byproduct in reaction SMARTS",
                            type=str2bool, default=False)
        parser.add_argument("--spectator", help="Add spectator in reaction SMARTS",
                            type=str2bool, default=False)
        parser.add_argument("--full", help="Get overall reaction instead of elementary steps",
                            type=str2bool, default=False)
        parser.add_argument("--end", help="Get termination reactions (products>>products)",
                            type=str2bool, default=True)
        parser.add_argument("--plain", help="Get reaction SMARTS without atom-mapping",
                            type=str2bool, default=False)
        parser.add_argument("--explicit_H", help="SMILES with Explicit H",
                            type=str2bool, default=True)
        parser.add_argument("--reagent",
                            help="Locate reagents at middle of the reaction (reactants>reagents>products) instead of both sides of reactants and products.",
                            type=str2bool, default=False)
        parser.add_argument("--remapping", help="Re-atom-mapping for each elementary reaction to start with 1.",
                            type=str2bool, default=True)

        # Arguments for data loading and saving
        parser.add_argument('--data', help='Path to the reaction data',
                            type=str, default='./data/test_data.txt')
        parser.add_argument('--save', help='Path to the saving data',
                            type=str, default='./results/test.txt')
        parser.add_argument("--debug", help="Collect all the failed reactions",
                            type=str2bool, default=False)
        parser.add_argument("--all_info",
                            help="Save all the information you can get, or you can get only elementary steps if not",
                            type=str2bool, default=False)
        parser.add_argument("--stat", help="Save the statistics dictionary",
                            type=str2bool, default=True)
        parser.add_argument("--rxn_class",
                            help="Specify the reaction class to work on. If not specified, it will work on all classes.",
                            type=str, default='')
        parser.add_argument("--process", help="The number of worker processes",
                            type=int, default=10)
        parser.add_argument("--verbosity",
                            help="control the verbosity; 0: silent, 1: prints the critical errors, 2: prints some details (do not recommend more than 10 reactions), 3: prints the process, 4: prints all",
                            type=int, default=0)

        return parser.parse_args()

    def get_mechanistic_network(rxn, args):
        conditions = calling_rxn_template(rxn)
        for i, cond in enumerate(conditions):
            reaction = reaction_process.Reaction_Network(rxn, i, args)
            reaction.set_template_dict(cond)
            for step in range(reaction.length):
                try:
                    reaction.run_reaction(step)
                except NoReactantError as e:
                    print(i, 'No reactant')
                    continue
                except NoAcidBaseError as e:
                    print(i, 'No acid base')
            if reaction.is_product_formed():
                elem_reactions = Get_Reactions(reaction)
                elem_reactions.print_graph()
                if not args.do_not_pruning:
                    elem_reactions.graph_pruning()
                elem_reactions.print_graph()
                elem_reactions.get_elementary_reactions_info()
                rxn_smi = elem_reactions.convert_to_smiles()
                print(rxn_smi)



    args = parse_arguments()
    # line = 'BrC1=C(Br)C=CC=C1C.OB(C)O.[Cl-].ClC1=CC=CC=C1>Cl[Pd]Cl.N.P(C2=CC=CC=C2)(C3=CC=CC=C3)C4=CC=CC=C4.P(C5=CC=CC=C5)(C6=CC=CC=C6)C7=CC=CC=C7.O>BrC8=C(C)C=CC=C8C Bromo Suzuki coupling'
    line = '[CH2:23]([Cl:24])[Cl:25].Cl.[NH:1]([C@:2]1([C:3]([NH:4][S:5]([CH:6]2[CH2:7][CH2:8]2)(=[O:9])=[O:10])=[O:11])[CH2:12][C@H:13]1[CH2:14][CH3:15])[C:18]([O:17][C:16]([CH3:20])([CH3:21])[CH3:22])=[O:19]>>[NH2:1][C@:2]1([C:3]([NH:4][S:5]([CH:6]2[CH2:7][CH2:8]2)(=[O:9])=[O:10])=[O:11])[CH2:12][C@H:13]1[CH2:14][CH3:15] N-Boc deprotection'
    rxn = line.split()[0]
    if len(line.split()[1:]) == 1:
        label = line.split()[-1]
    else:
        label = ' '.join(line.split()[1:])

    rxn_dict = {'reaction_name': label,
                'reaction_smiles': rxn, }
    rxn_dict = reagent_matching_for_single_reaction(rxn_dict, label)

    get_mechanistic_network(rxn_dict, args)
