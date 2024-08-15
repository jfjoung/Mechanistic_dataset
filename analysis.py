import os
import pickle
from scripts.analysis_setting import Args
from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import defaultdict
from templates import AcidBase_lookup, Reaction_templates
from multiprocessing import Pool
from tqdm import tqdm
import pandas as pd
from utils.validity_check import check_reaction_validity

def loadall(filename):
    with open(filename, "rb") as f:
        while True:
            try:
                yield pickle.load(f)
            except EOFError:
                break

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
        symbol1, symbol2 = mapping_to_type[int(atom1)], mapping_to_type[int(atom2)]
        sorted_pair = sorted([symbol1, symbol2], key=lambda x: (priority.get(x, 7), x))
        bond_change_types[''.join(sorted_pair)] += 1

    for atom1, _ in hydrogen_changes:
        symbol1 = mapping_to_type[int(atom1)]
        sorted_pair = sorted([symbol1, 'H'], key=lambda x: (priority.get(x, 7), x))
        bond_change_types[''.join(sorted_pair)] += 1

    bond_edits = bond_changes.union(hydrogen_changes)
    
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

    return ers, bond_change_types

def merge_dicts(d, u):
    for k, v in u.items():
        if isinstance(v, dict):
            d[k] = merge_dicts(d.get(k, {}), v)
        else:
            d[k] = d.get(k, 0) + v
    return d

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
    print(f"Covered acid-base pair: {len(AB)}")

def reaction_analysis(input):

    reaction, args = input 


    # Dictionaries
    # ERS
    ers_dict = {}
    bond_change_types = {}
    rp_dict = {}

    #Mol
    smiles_identity = {}
    MW_dict = {}
    Atom_dict = {}
    Atom_types = {}

    invalid_reaction = None

    
    reaction_info = reaction.reaction_info

    if args.ERS_analysis:
        for pathways in reaction_info.values():
            for rxn_i, chemicals in pathways.items():
                rnodes = chemicals['reactants']
                pnodes = chemicals['products']
                if rnodes != pnodes:
                    try:
                        rp_count = f"{len(rnodes)}r{len(pnodes)}p"
                        if rp_count in rp_dict:
                            rp_dict[rp_count] += 1
                        else:
                            rp_dict[rp_count] = 1

                        rsmi = '.'.join([reaction.get_node(node_id).smiles_w_mapping for node_id in rnodes])
                        psmi = '.'.join([reaction.get_node(node_id).smiles_w_mapping for node_id in pnodes])
                        rxn_smi = '>>'.join([rsmi, psmi])
                        ers, bc = get_bond_edits(rxn_smi)
                        merge_dicts(bond_change_types, bc)

                        if ers in ers_dict:
                            ers_dict[ers] += 1
                        else:
                            ers_dict[ers] = 1
                    except:
                        continue

    if args.mol_analysis:
        chemical_nodes = reaction.reactant_node + reaction.product_node + reaction.byproduct_node + reaction.intermediate_node + reaction.spectator_node

        for node_id in chemical_nodes:
            smi = reaction.get_node(node_id).smiles
            identity = reaction.get_node(node_id).identity
            if smi not in smiles_identity:
                smiles_identity[smi] = {identity: 1}
            else:
                if identity in smiles_identity[smi]:
                    smiles_identity[smi][identity] += 1
                else:
                    smiles_identity[smi][identity] = 1

            mol = Chem.MolFromSmiles(smi, sanitize=False)
            mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(mol,
                            Chem.SanitizeFlags.SANITIZE_FINDRADICALS | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                            catchErrors=True)
            
            mass = round(Descriptors.MolWt(mol))
            atom_count = mol.GetNumHeavyAtoms()

            atom_type = {}
            for a in mol.GetAtoms():
                symbol = a.GetSymbol()
                if symbol not in atom_type:
                    atom_type[symbol] = 1
                else:
                    atom_type[symbol] += 1
            merge_dicts(Atom_types, atom_type)


            mw_dict = {identity: {mass: 1}}
            atom_dict = {identity: {atom_count: 1}}
            merge_dicts(MW_dict, mw_dict)
            merge_dicts(Atom_dict, atom_dict)

    if args.validity:
        for rxn_smi in reaction.rxn_smi:
            # print(rxn_smi)
            invalid = False
            try:
                result = check_reaction_validity(rxn_smi)
                if not result:
                    if args.verbosity:
                        reaction.print_graph()
                        print(rxn_smi)
                        print(reaction_info)
                    invalid = True
            except Exception as e:
                if args.verbosity:
                    reaction.print_graph()
                    print(rxn_smi)
                    print(e)
                    print(reaction_info)
                invalid = True
            if invalid:
                reaction_smiles = reaction.reaction_smiles
                reaction_class = reaction.reaction_class
                invalid_reaction = f'{reaction_smiles} {reaction_class}'
                break

    return ers_dict, rp_dict, bond_change_types, smiles_identity, MW_dict, Atom_dict, Atom_types, invalid_reaction

        

def main(args):
    save_file = args.data

    base_file_root, _ = os.path.splitext(args.data)
    ers_file_path = f"{base_file_root}_ERS.csv"
    rp_file_path = f"{base_file_root}_RnPn.csv"
    bc_file_path = f"{base_file_root}_BondChange.csv"

    MW_file_path = f"{base_file_root}_MW.csv"
    AtomCount_file_path = f"{base_file_root}_AtomCount.csv"
    AtomTypes_file_path = f"{base_file_root}_AtomTypes.csv"

    validity_file_path = f"{base_file_root}_invalid.txt"
   

    if args.template_analysis:
        template_analysis(args)

    items = loadall(save_file)
    ers_dict = {}
    bond_change_types = {}
    rp_dict = {}

    smiles_identity = {}
    MW_dict = {}
    Atom_dict = {}
    Atom_types = {}

    invalid_reactions = []
    p = Pool(args.process)

    iterables = [(rxn, args) for reactions in tqdm(items, desc="Data loading") for rxn in reactions]

    for results in tqdm(p.imap_unordered(reaction_analysis, iterables), total=len(iterables)):
        ers, rp, bc, iden, mw, atom_count, atom_type, invalid = results

        merge_dicts(ers_dict, ers)
        merge_dicts(rp_dict, rp)
        merge_dicts(bond_change_types, bc)
        merge_dicts(smiles_identity, iden)
        merge_dicts(MW_dict, mw)
        merge_dicts(Atom_dict, atom_count)
        merge_dicts(Atom_types, atom_type)
        if invalid is not None:
            invalid_reactions.append(invalid)

    if args.ERS_analysis:
        ers_df = pd.DataFrame(list(ers_dict.items()), columns=['Types', 'Count']).sort_values(by='Types')
        ers_df.to_csv(ers_file_path, index=False)

        rp_df = pd.DataFrame(list(rp_dict.items()), columns=['Types', 'Count']).sort_values(by='Types')
        rp_df.to_csv(rp_file_path, index=False)

        bc_df = pd.DataFrame(list(bond_change_types.items()), columns=['Types', 'Count']).sort_values(by='Count', ascending=False)
        bc_df.to_csv(bc_file_path, index=False)

    if args.mol_analysis:
        print(f'There are {len(smiles_identity)} unique molecules')

        MW_df = pd.DataFrame.from_dict(MW_dict, orient='index')
        MW_df = MW_df.transpose().fillna(0).sort_index()
        MW_df.to_csv(MW_file_path, index_label="Molecular weight")

        atom_df = pd.DataFrame.from_dict(Atom_dict, orient='index')
        atom_df = atom_df.transpose().fillna(0).sort_index()
        atom_df.to_csv(AtomCount_file_path, index_label="Number of heavy atoms")

        atomtype_df =pd.DataFrame(list(Atom_types.items()), columns=['Types', 'Count']).sort_values(by='Count', ascending=False)
        atomtype_df.to_csv(AtomTypes_file_path, index=False)

        print(f'There are {len(Atom_types)} types of atoms')

    if args.validity:
        print(f'There are {len(invalid_reactions)} invalid reactions')
        with open(validity_file_path, 'w') as fout:
            fout.write('\n'.join(invalid_reactions) + '\n')



if __name__ == '__main__':
    args=Args()
    main(args)