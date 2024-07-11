import re

def parse_smarts(smarts):
    # Extended pattern to capture complex atoms with additional properties like ;H1, $([O]-S), etc.
    # Modified to include '*' as a valid atom symbol.
    atom_pattern = re.compile(r'\[([*#A-Za-z0-9,+;\-&!$()\[\]]*):(\d+)\]')
    atoms = atom_pattern.findall(smarts)
    atom_dict = {}

    for atom, mapping in atoms:
        if 'H' in atom:
            h_count = re.search(r'H(\d*)', atom)
            if h_count and h_count.group(1):
                h_count = int(h_count.group(1))
            else:
                h_count = 1
        else:
            h_count = 0

        atom_dict[mapping] = {'atom': atom, 'h_count': h_count}

    return atom_dict

def clean_hydrogen_notation(atom):
    atom = re.sub(r'H(\d+)', '', atom)
    atom = re.sub(r'H(?!0)', '', atom)
    atom = re.sub(r';{2,}', ';', atom).replace(';:', ':').replace(';,', ';')
    return atom

def update_smarts_with_hydrogen_mapping(smarts, atom_dict):
    used_mappings = set(int(key) for key in atom_dict.keys())
    new_smarts = smarts
    added_h_mappings = {}

    for mapping, info in atom_dict.items():
        atom_str = f'[{info["atom"]}:{mapping}]'
        cleaned_atom_str = clean_hydrogen_notation(atom_str)

        new_h_atoms = []
        for _ in range(info['h_count']):
            new_mapping = min(set(range(1, len(atom_dict) * 3)) - used_mappings)
            used_mappings.add(new_mapping)
            new_h_atoms.append(f'([H:{new_mapping}])')
            if mapping not in added_h_mappings:
                added_h_mappings[mapping] = []
            added_h_mappings[mapping].append(new_mapping)

        new_atom_str = cleaned_atom_str + ''.join(new_h_atoms)
        new_smarts = new_smarts.replace(atom_str, new_atom_str)

    for mapping in atom_dict.keys():
        if mapping not in added_h_mappings:
            added_h_mappings[mapping] = None

    return new_smarts, added_h_mappings

def detect_hydrogen_changes(reactant_dict, product_dict):
    hydrogen_changes = {'lost': {}, 'gained': {}}
    total_lost = 0
    total_gained = 0

    for key in reactant_dict:
        if key in product_dict:
            r_h_count = reactant_dict[key]['h_count']
            p_h_count = product_dict[key]['h_count']
            if r_h_count > p_h_count:
                hydrogen_changes['lost'][key] = r_h_count - p_h_count
                total_lost += r_h_count - p_h_count
            elif r_h_count < p_h_count:
                hydrogen_changes['gained'][key] = p_h_count - r_h_count
                total_gained += p_h_count - r_h_count

    return hydrogen_changes, total_lost, total_gained

def update_smarts_product(smarts, atom_dict, added_h_mappings, hydrogen_changes):
    used_mappings = set(int(key) for key in atom_dict.keys())
    for mapping_list in added_h_mappings.values():
        if mapping_list is not None:
            used_mappings.update(mapping_list)

    new_smarts = smarts

    for mapping, info in atom_dict.items():
        atom_str = f'[{info["atom"]}:{mapping}]'
        cleaned_atom_str = clean_hydrogen_notation(atom_str)

        new_h_atoms = []
        if mapping in hydrogen_changes['gained']:
            gain_count = hydrogen_changes['gained'][mapping]
            for _ in range(gain_count):
                if hydrogen_changes['lost']:
                    lost_mapping = next(iter(hydrogen_changes['lost']))
                    if added_h_mappings[lost_mapping]:
                        new_mapping = added_h_mappings[lost_mapping].pop(0)
                        new_h_atoms.append(f'([H:{new_mapping}])')
                        hydrogen_changes['lost'][lost_mapping] -= 1
                        if hydrogen_changes['lost'][lost_mapping] == 0:
                            del hydrogen_changes['lost'][lost_mapping]

        remaining_h_count = info['h_count']
        if mapping in hydrogen_changes['gained']:
            remaining_h_count -= hydrogen_changes['gained'][mapping]
        for _ in range(remaining_h_count):
            if added_h_mappings.get(mapping):
                new_mapping = added_h_mappings[mapping].pop(0)
                new_h_atoms.append(f'([H:{new_mapping}])')

        new_atom_str = cleaned_atom_str + ''.join(new_h_atoms)
        new_smarts = new_smarts.replace(atom_str, new_atom_str)

    return new_smarts

def update_smarts_product_no_change(smarts, atom_dict, added_h_mappings):
    new_smarts = smarts

    for mapping, info in atom_dict.items():
        atom_str = f'[{info["atom"]}:{mapping}]'
        cleaned_atom_str = clean_hydrogen_notation(atom_str)

        new_h_atoms = []
        if mapping not in added_h_mappings:
            added_h_mappings[mapping] = []

        if added_h_mappings[mapping]:
            for h_mapping in added_h_mappings[mapping]:
                new_h_atoms.append(f'([H:{h_mapping}])')

        new_atom_str = cleaned_atom_str + ''.join(new_h_atoms)
        new_smarts = new_smarts.replace(atom_str, new_atom_str)

    return new_smarts

def modify_explicit_H(smarts):
    reactants, reagent, products = smarts.split('>')
    reactant_dict = parse_smarts(reactants)
    product_dict = parse_smarts(products)

    if all(reactant_dict[key]['h_count'] == product_dict[key]['h_count'] for key in reactant_dict if key in product_dict):
        new_reactants, added_h_mappings = update_smarts_with_hydrogen_mapping(reactants, reactant_dict)
        new_products = update_smarts_product_no_change(products, product_dict, added_h_mappings)
    else:
        hydrogen_changes, total_lost, total_gained = detect_hydrogen_changes(reactant_dict, product_dict)
        if total_lost != total_gained:
            raise ValueError(f"Total lost hydrogens ({total_lost}) does not match total gained hydrogens ({total_gained}).")
        new_reactants, added_h_mappings = update_smarts_with_hydrogen_mapping(reactants, reactant_dict)
        new_products = update_smarts_product(products, product_dict, added_h_mappings, hydrogen_changes)

    return f'{new_reactants}>{reagent}>{new_products}', added_h_mappings


if __name__ == '__main__':
    # Test case 1: No hydrogen change
    smarts1 = '[Cl,Br,I,O:2]-[Pd;+0:1]-[P;H0;+1:4].[OH2;+0:5]>>[Cl,Br,I,O:2]-[Pd;+0:1]-[P;H0;+0:4]-[OH2;+1:5]'
    result1, added_h_mappings1 = modify_explicit_H(smarts1)
    print(result1)

    # Test case 2: Hydrogen change
    smarts2 = '[#7;H3;+0:1].[#8;H0:8]-[Pd:6]-[#6;+0:5]>>[#7;H2;+0:1]-[Pd:6]-[#6;+0:5].[#8;H1;+0:8]'
    result2, added_h_mappings2 = modify_explicit_H(smarts2)
    print(result2)

    # Test case 3
    smarts3 = '[Cl,Br,I,O:2]-[Pd;H1;+0:1]>>[Cl,Br,I,O;H1:2].[Pd;H0;+0:1]'
    result3, added_h_mappings3 = modify_explicit_H(smarts3)
    print(result3)

    # Test case 4
    smarts4 = '[#8;H2:8].[Cl,Br,I,O$([O]-S):7]-[Pd:6]-[#6;+0:5]>>[Cl,Br,I,O$([O]-S);H1;+0:7].[#8;H1;+0:8]-[Pd:6]-[#6;+0:5]'
    result4, added_h_mappings4 = modify_explicit_H(smarts4)
    print(result4)

    # Test case 5: No hydrogen
    smarts5 = '[Cl,Br,I,O$([O]-S):7]-[#6;+0:5].[Pd&!$([Pd]-I)&!$([Pd]-Br)&!$([Pd]-Cl):6]>>[Cl,Br,I,O$([O]-S):7]-[Pd:6]-[#6;+0:5]'
    result5, added_h_mappings5 = modify_explicit_H(smarts5)
    print(result5)

    # Test cast 6: Hydrogen on aromatic ring
    smarts6 = '[*:1][O:2][N+:3]([O-:4])=[O:5].[*;a;+0:8]:[c;H1;+0:9]>>[*;a;+1:8]:[c;H1;+0:9]-[N+:3](-[O-:4])(=[O:5]).[*:1][O;H0;-1:2]'
    result6, added_h_mappings6 = modify_explicit_H(smarts6)
    print(result6)

    # Test case 7
    smarts7 = '[O;H1;+1:1]=[C;H0;+0:2]-[O;H0;+0:4]-[#6;+0:5].[n&!$([#7]~[#7]);H0:100]1[c:99][c:98][c:97][c:96][c:95]1>>[O;H0;+0:1]=[C;H0;+0:2]-[O;H0;+0:4]-[#6;+0:5].[n;+1&!$([#7]~[#7]);H1:100]1[c:99][c:98][c:97][c:96][c:95]1'
    result7, added_h_mappings7 = modify_explicit_H(smarts7)
    print(result7)