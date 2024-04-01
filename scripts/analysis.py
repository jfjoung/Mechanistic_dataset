from multiprocessing import Pool
import jsonlines
import json
import networkx as nx
from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import defaultdict
import matplotlib.pyplot as plt
from tqdm import tqdm
import os
from utils.apply_mechanistic_template import find_chemical_nodes


def mass_analysis(mol):
    mass = Descriptors.MolWt(mol)
    return mass

def heavy_atom_count(mol):
    return mol.GetNumHeavyAtoms()


def merge_dicts(d, u):
    for k, v in u.items():
        if isinstance(v, dict):
            d[k] = merge_dicts(d.get(k, {}), v)
        else:
            d[k] = d.get(k, 0) + v
    return d

def node_analysis(G, nid, MW_dict, Atom_dict, args):
    smi = G.nodes[nid]['molecule']['smiles_w_mapping']
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol,
                     Chem.SanitizeFlags.SANITIZE_FINDRADICALS | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                     catchErrors=True)

    molwt = mass_analysis(mol)
    bin_key = (int(molwt) // args.weight_bin) * args.weight_bin
    MW_dict[bin_key] += 1

    num_atom = heavy_atom_count(mol)
    bin_key = (int(num_atom) // args.atom_bin) * args.atom_bin
    Atom_dict[bin_key] += 1

def analysis_single(input):
    rxn_dict, args = input

    reactant_MW_dict = defaultdict(int)
    reactant_Atom_dict = defaultdict(int)

    product_MW_dict = defaultdict(int)
    product_Atom_dict = defaultdict(int)

    inter_MW_dict = defaultdict(int)
    inter_Atom_dict = defaultdict(int)

    by_MW_dict = defaultdict(int)
    by_Atom_dict = defaultdict(int)

    spec_MW_dict = defaultdict(int)
    spec_Atom_dict = defaultdict(int)

    rxn_network = rxn_dict['Mechanism']
    for rxn_network in rxn_network.values():
        G = nx.node_link_graph(rxn_network['Reaction graph'])
        reaction_nodes, reactant_nodes, product_nodes, byproduct_nodes, intermediate_nodes, spectator_nodes = find_chemical_nodes(G)

        for nid in reactant_nodes:
            node_analysis(G, nid, reactant_MW_dict, reactant_Atom_dict, args)
        for nid in product_nodes:
            node_analysis(G, nid, product_MW_dict, product_Atom_dict, args)
        for nid in intermediate_nodes:
            node_analysis(G, nid, inter_MW_dict, inter_Atom_dict, args)
        for nid in byproduct_nodes:
            node_analysis(G, nid, by_MW_dict, by_Atom_dict, args)
        for nid in spectator_nodes:
            node_analysis(G, nid, spec_MW_dict, spec_Atom_dict, args)

    molecular_weight_dict = {'Reactant': dict(sorted(reactant_MW_dict.items())),
                             'Product': dict(sorted(product_MW_dict.items())),
                             'Intermediate': dict(sorted(inter_MW_dict.items())),
                             'Byproduct': dict(sorted(by_MW_dict.items())),
                             'Spectator': dict(sorted(spec_MW_dict.items())),
                             }
    heavy_atom_dict = {'Reactant': dict(sorted(reactant_Atom_dict.items())),
                       'Product': dict(sorted(product_Atom_dict.items())),
                       'Intermediate': dict(sorted(inter_Atom_dict.items())),
                       'Byproduct': dict(sorted(by_Atom_dict.items())),
                       'Spectator': dict(sorted(spec_Atom_dict.items())),
                       }

    return molecular_weight_dict , heavy_atom_dict

def drawing_fig(distributions, args, file_name='MW'):
    # Initialize figure and axes for the bar graph
    fig, ax = plt.subplots(figsize=(12, 8))
    # Assign a unique color to each chemical species
    colors = ['blue', 'orange', 'green', 'red', 'purple']

    # # Collect all unique keys (MW ranges) from all categories to create a comprehensive set of x-axis labels
    # all_keys = set(key for dist in distribution.values() for key in dist)
    # sorted_keys = sorted(all_keys)
    #
    # # Plotting each chemical species separately
    # for i, (species, counts) in enumerate(distribution.items()):
    #     # Adjusting positions for separation
    #     positions = [x + (i - len(distribution) / 2) * 0.1 for x in range(len(sorted_keys))]
    #
    #     # Only plotting bars for keys present in this species
    #     heights = [counts[k] if k in counts else 0 for k in sorted_keys]
    #
    #     ax.bar(positions, heights, width=0.1, label=species, alpha=0.75, color=colors[i % len(colors)])
    #
    # # Setting the x-axis ticks to show the correct MW range
    # ax.set_xticks(range(len(sorted_keys)))
    #
    # if file_name == 'MW':
    #     ax.set_xticklabels([f"{k}-{k + args.weight_bin - 1}" for k in sorted_keys], rotation=45)
    # elif file_name == 'Atom':
    #     ax.set_xticklabels([f"{k}-{k + args.atom_bin - 1}" for k in sorted_keys], rotation=45)
    # Define the full range of x values (number of heavy atoms) across all categories for completeness

    all_keys = set().union(*[distribution.keys() for distribution in distributions.values()])
    full_range_keys = sorted(all_keys)

    for i, (category, distribution) in enumerate(distributions.items()):
        # Filling missing keys with 0 values for a complete distribution across the full range
        complete_values = [distribution.get(key, 0) for key in full_range_keys]

        ax.plot(full_range_keys, complete_values, label=category, linestyle='-', color=colors[i])
        ax.fill_between(full_range_keys, 0, complete_values, color=colors[i], alpha=0.3)

    # Adding labels and title
    ax.set_xlabel(f'{file_name} Range')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{file_name} Distribution by Chemical Species')
    if args.log_scale:
        ax.set_yscale('log')
    plt.legend()

    # Saving the figure
    plt.tight_layout()

    base_file_root, _ = os.path.splitext(args.data)
    stat_file_path = f"{base_file_root}_{file_name}.png"
    plt.savefig(stat_file_path)

def analysis(args):
    fpath = args.data

    p = Pool(args.process)

    molecular_weight_dict = {}
    heavy_atom_dict = {}
    with jsonlines.open(fpath) as f:
        iterables = [(rxn_dict, args) for rxn_dict in f.iter()]
        for result in tqdm(p.imap(analysis_single, iterables), total=len(iterables)):
            MW_dict, Atom_dict = result
            merge_dicts(molecular_weight_dict, MW_dict)
            merge_dicts(heavy_atom_dict, Atom_dict)


    molecular_weight_dict = {key: dict(sorted(value.items())) for key, value in molecular_weight_dict.items()}
    heavy_atom_dict = {key: dict(sorted(value.items())) for key, value in heavy_atom_dict.items()}
    print(f'Molecular weight distribution = {molecular_weight_dict}')
    print(f'The number of heavy atom distribution = {heavy_atom_dict}')

    drawing_fig(molecular_weight_dict, args, file_name='MW')
    drawing_fig(heavy_atom_dict, args, file_name='Atom')

def molecule_single(input):
    rxn_dict, args = input
    rxn_network = rxn_dict['Mechanism']
    smiles_identity_count_updated = {}
    for rxn_network in rxn_network.values():
        G = nx.node_link_graph(rxn_network['Reaction graph'])
        reaction_nodes, reactant_nodes, product_nodes, byproduct_nodes, intermediate_nodes, spectator_nodes = find_chemical_nodes(G)
        chemical_nodes = reactant_nodes+product_nodes+byproduct_nodes+intermediate_nodes+spectator_nodes

        for nid in chemical_nodes:
            item = G.nodes[nid]

            smiles = item['molecule']['smiles']
            identity = item['molecule']['identity']

            if smiles not in smiles_identity_count_updated:
                smiles_identity_count_updated[smiles] = {identity: 1}
            else:
                if identity in smiles_identity_count_updated[smiles]:
                    smiles_identity_count_updated[smiles][identity] += 1
                else:
                    smiles_identity_count_updated[smiles][identity] = 1
    return smiles_identity_count_updated

def mol_analysis(args):
    fpath = args.data

    p = Pool(args.process)

    molecule_dict = {}
    with jsonlines.open(fpath) as f:
        iterables = [(rxn_dict, args) for rxn_dict in f.iter()]
        for result in tqdm(p.imap(molecule_single, iterables), total=len(iterables)):
            merge_dicts(molecule_dict, result)

    print(f'There are {len(molecule_dict)} unique molecules')

    base_file_root, _ = os.path.splitext(args.data)
    mol_file_path = f"{base_file_root}_molecule_dict.json"
    with open(mol_file_path, 'w') as file:
        json.dump(molecule_dict, file, indent=4)





