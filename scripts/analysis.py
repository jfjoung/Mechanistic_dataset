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
from utils import AcidBase_lookup, Reaction_templates


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
        for result in tqdm(p.imap_unordered(analysis_single, iterables), total=len(iterables)):
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
        for result in tqdm(p.imap_unordered(molecule_single, iterables), total=len(iterables)):
            merge_dicts(molecule_dict, result)

    print(f'There are {len(molecule_dict)} unique molecules')

    base_file_root, _ = os.path.splitext(args.data)
    mol_file_path = f"{base_file_root}_molecule_dict.json"
    with open(mol_file_path, 'w') as file:
        json.dump(molecule_dict, file, indent=4)

def get_bond_edits(rxn_smi):
    reactants = Chem.MolFromSmiles(rxn_smi.split('>')[0], sanitize=False)
    products = Chem.MolFromSmiles(rxn_smi.split('>')[2], sanitize=False)

    conserved_maps = [a.GetIntProp('molAtomMapNumber') for a in reactants.GetAtoms() if
                      a.GetIntProp('molAtomMapNumber')]
    conserved_maps_p = [a.GetIntProp('molAtomMapNumber') for a in products.GetAtoms() if
                      a.GetIntProp('molAtomMapNumber')]

    mapping_to_type = {}
    for a in reactants.GetAtoms():
        mapping_to_type[a.GetAtomMapNum()] = a.GetSymbol()

    bond_changes = set()
    conserved_maps_a = [a.GetAtomMapNum() for a in products.GetAtoms() if a.HasProp('molAtomMapNumber')]
    hydrogen_changes = set()
    bonds_prev = {}
    for bond in reactants.GetBonds():
        nums = sorted(
            [bond.GetBeginAtom().GetIntProp('molAtomMapNumber'), bond.GetEndAtom().GetIntProp('molAtomMapNumber')])
        if (nums[0] not in conserved_maps) or (nums[1] not in conserved_maps_p): continue
        bonds_prev['{}~{}'.format(nums[0], nums[1])] = bond.GetBondTypeAsDouble()

    bonds_new = {}
    for bond in products.GetBonds():
        nums = sorted(
            [bond.GetBeginAtom().GetIntProp('molAtomMapNumber'), bond.GetEndAtom().GetIntProp('molAtomMapNumber')])
        if (nums[0] not in conserved_maps) or (nums[1] not in conserved_maps): continue
        bonds_new['{}~{}'.format(nums[0], nums[1])] = bond.GetBondTypeAsDouble()

    hydrogen_prev = {}
    for atom in reactants.GetAtoms():
        if atom.GetAtomMapNum() not in conserved_maps_a: continue
        hydrogen_prev['{}'.format(atom.GetAtomMapNum())] = atom.GetTotalNumHs()

    hydrogen_new = {}
    for atom in products.GetAtoms():
        hydrogen_new['{}'.format(atom.GetAtomMapNum())] = atom.GetTotalNumHs()

    for bond in bonds_prev:
        if bond not in bonds_new:
            bond_changes.add((bond.split('~')[0], bond.split('~')[1], -bonds_prev[bond]))  # lost bond
        else:
            if bonds_prev[bond] != bonds_new[bond]:
                bond_changes.add(
                    (bond.split('~')[0], bond.split('~')[1], bonds_new[bond] - bonds_prev[bond]))  # changed bond

    for bond in bonds_new:
        if bond not in bonds_prev:
            bond_changes.add((bond.split('~')[0], bond.split('~')[1], bonds_new[bond]))  # new bond

    for atom in hydrogen_new:
        if hydrogen_prev[atom] != hydrogen_new[atom]:
            hydrogen_changes.add((atom, hydrogen_new[atom] - hydrogen_prev[atom]))  # changed hydrogen

    bond_change_types = defaultdict(int)
    priority = {'C': 1, 'N': 2, 'O': 3, 'P': 4, 'S': 5, 'H': 6}

    for atom1, atom2, _ in bond_changes:
        # 원자 매핑 번호를 기반으로 원자 기호 찾기
        symbol1, symbol2 = mapping_to_type[int(atom1)], mapping_to_type[int(atom2)]

        # 원자 기호 쌍 정렬 (CNOPS 순서, H 맨 뒤)
        sorted_pair = sorted([symbol1, symbol2], key=lambda x: (priority.get(x, 7), x))

        # 정렬된 쌍을 문자열로 결합하여 카운트
        bond_change_types[''.join(sorted_pair)] += 1

    for atom1, _ in hydrogen_changes:
        symbol1 = mapping_to_type[int(atom1)]
        sorted_pair = sorted([symbol1, 'H'], key=lambda x: (priority.get(x, 7), x))
        bond_change_types[''.join(sorted_pair)] += 1


    return bond_changes.union(hydrogen_changes), bond_change_types

def get_elementary_reaction_step(input):
    rxn_dict, args = input
    ers_dict = {}
    bond_change_types = {}
    rp_dict = {}
    for key in rxn_dict['Mechanism'].keys():
        rxns = rxn_dict['Mechanism'][key]['Elementary steps']



        for rxn in rxns:
            try:
                bond_edits, bond_change_type = get_bond_edits(rxn)
                merge_dicts(bond_change_types, bond_change_type)

                r_smi, p_smi = rxn.split('>>')
                if '.' in r_smi:
                    r_count = 2
                else:
                    r_count = 1

                p_count = len(p_smi.split('.'))

                rp_count = f"{r_count}r{p_count}p"

                if rp_count in rp_dict:
                    rp_dict[rp_count] += 1
                else:
                    rp_dict[rp_count] = 1

            except: continue
            positive_count, negative_count = 0, 0

            for item in bond_edits:
                if item[-1] >= 1:
                    positive_count += 1
                elif item[-1] <= -1:
                    negative_count += 1

            ers = ""
            if negative_count > 0:
                ers += f"b{negative_count}"
            if positive_count > 0:
                ers += f"f{positive_count}"

            # if ers == 'b2f3':
            #     print(rxn)
            #     print(bond_edits)
            #     print(ers)

            if ers in ers_dict:
                ers_dict[ers] += 1
            else:
                ers_dict[ers] = 1

    return ers_dict, bond_change_types, rp_dict


def ERS_analysis(args):
    fpath = args.data
    p = Pool(args.process)

    ers_dict = {}
    bond_change_types = {}
    rp_dict = {}
    with jsonlines.open(fpath) as f:
        iterables = [(rxn_dict, args) for rxn_dict in f.iter()]
        for result in tqdm(p.imap_unordered(get_elementary_reaction_step, iterables), total=len(iterables)):
        # for input in iterables:
        #     out_dict = get_elementary_reaction_step(input)
            ers, bct, rp = result
            merge_dicts(ers_dict, ers)
            merge_dicts(bond_change_types, bct)
            merge_dicts(rp_dict, rp)

    def sort_key(item):
        key, _ = item
        # b와 f의 존재 여부를 판단하여 정렬 우선순위 결정
        b_exists = 'b' in key
        f_exists = 'f' in key
        if b_exists and not f_exists:
            # b만 존재하는 경우
            return (1, key)
        elif not b_exists and f_exists:
            # f만 존재하는 경우
            return (2, key)
        else:
            # b와 f 둘 다 존재하는 경우
            return (3, key)

    ers_dict = dict(sorted(ers_dict.items(), key=sort_key))
    print(ers_dict)
    bond_change_types = dict(sorted(bond_change_types.items(), key=lambda item: (-item[1], item[0])))
    print(bond_change_types)
    print(rp_dict)


    a = 1.5

    plt.rcParams.update({'font.size': 5 * a})
    fig = plt.figure(figsize=(5 * a, 3 * a), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[13, 10], height_ratios=[1, 1])

    labels = list(bond_change_types.keys())[:20]
    values = list(bond_change_types.values())[:20]
    ax1 = fig.add_subplot(gs[0, :])
    bars = ax1.bar(labels, values, color='skyblue', width=0.3)
    ax1.set_yscale('log')
    ax1.tick_params(labelrotation=45)

    digits_min = len(str(min(values)))-1
    digits_max = len(str(max(values)))
    ax1.set_ylim([10**digits_min, 3*10**digits_max])

    ax1.set_xlabel('Bond change type')  # , fontsize=14)
    ax1.set_ylabel('$\\mathbf{N}_{reaction}$')  # , fontsize=14)
    ax1.text(0.95, 0.95, 'Total 81 types of bond changes', horizontalalignment='right', verticalalignment='top',
             transform=ax1.transAxes)

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2.0, height, f'{height}', ha='center', va='bottom', rotation=20)

    labels = list(ers_dict.keys())
    values = list(ers_dict.values())

    ax_bottom_left = fig.add_subplot(gs[1, 0])
    bars = ax_bottom_left.bar(labels, values, color='skyblue', width=0.3)
    ax_bottom_left.set_yscale('log')
    ax_bottom_left.tick_params(labelrotation=45)
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2.0, height, f'{height}', ha='center', va='bottom', rotation=0)

    digits_min = len(str(min(values)))-1
    digits_max = len(str(max(values)))
    ax_bottom_left.set_ylim([10**digits_min, 3*10**digits_max])

    # ax_bottom_left.set_ylim([50, 3000000])
    ax_bottom_left.set_xlabel('Reaction type')  # , fontsize=14)
    ax_bottom_left.set_ylabel('$\\mathbf{N}_{reaction}$')  # , fontsize=14)
    ax_bottom_left.set_title(' ')

    labels = list(rp_dict.keys())
    values = list(rp_dict.values())
    ax_bottom_right = fig.add_subplot(gs[1, 1])
    bars = ax_bottom_right.bar(labels, values, color='skyblue', width=0.3)
    ax_bottom_right.set_yscale('log')
    ax_bottom_right.tick_params(labelrotation=45)

    digits_min = len(str(min(values)))-1
    digits_max = len(str(max(values)))
    ax_bottom_right.set_ylim([10**digits_min, 3*10**digits_max])

    # ax_bottom_right.set_ylim([10, 2000000])
    ax_bottom_right.set_xlim([-0.9, 9.9])

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2.0, height, f'{height}', ha='center', va='bottom', rotation=0)
    ax_bottom_right.set_xlabel('Reaction type')  # , fontsize=14)
    ax_bottom_right.set_ylabel('$\\mathbf{N}_{reaction}$')  # , fontsize=14)
    ax_bottom_right.set_title(' ')

    # Saving the figure
    plt.tight_layout()

    base_file_root, _ = os.path.splitext(args.data)
    stat_file_path = f"{base_file_root}_reaction_analysis.png"
    plt.savefig(stat_file_path)

def template_analysis(args):
    AB = AcidBase_lookup.Acid_base

    reaction_name = []
    for rxnclass in Reaction_templates.class_reaction_templates.keys():
        for rxn in rxnclass:
            reaction_name.append(rxn)

    print(f"Covered reaction classes: {len(reaction_name)}")

    templates = []
    condition = 0

    for rxnclass, value in Reaction_templates.class_reaction_templates.items():
        for rxntype in value.values():
            condition += 1
            for elemrxn in rxntype['Stages'].values():
                if elemrxn['Templates']:
                    for template in elemrxn['Templates']:
                        templates.append(template)

    print(f"Covered reaction conditions: {condition}")
    print(f"Unique reaction templates: {len(set(templates))}")