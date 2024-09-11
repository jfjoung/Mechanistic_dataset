from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

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

def atommap_to_isotope(mol):
    mol_copy = Chem.Mol(mol)
    for a in mol_copy.GetAtoms(): #Change atom mapping to isotope labeling
        a.SetIsotope(a.GetAtomMapNum())
        a.SetAtomMapNum(0)
    return mol_copy

class Molecule_Node:
    '''
    Process the molecules in the reaction flask
    '''
    def __init__(self, mol, args):
        self.mol = mol  #It has isotope labeling
        self.smiles_w_isotope = None
        self.atom_mapping = None
        self.smiles_w_mapping = None
        self.smiles = None
        self.identity = None
        self.args = args
        self.index = None
        self.check_validity()
        self.get_plain_smiles()

    def __str__(self):
        return self.smiles

    def check_validity(self):
        '''
        Check the chemical validity.
        Valance of the atom is not checked, because some intermediate might violate the octet rule.
        '''
        _mol = Chem.Mol(self.mol)
        try:
            _mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(_mol,
                             Chem.SanitizeFlags.SANITIZE_FINDRADICALS | Chem.SanitizeFlags.SANITIZE_SETAROMATICITY | Chem.SanitizeFlags.SANITIZE_SETCONJUGATION | Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION | Chem.SanitizeFlags.SANITIZE_SYMMRINGS,
                             catchErrors=True)
        except:
            raise ValueError("Not valid molecule")

    def get_plain_smiles(self):
        '''
        Convert mol object to plain SMILES, no atom mapping and no isotope labeling
        '''
        _mol = Chem.Mol(self.mol)
        plain = Chem.MolToSmiles(remove_atom_map(_mol, isotope=True))
        self.smiles = plain

    def add_reactant(self, idx):
        '''
        It gives two smiles labeled with atom mapping or isotope.
        It is for the reactants.
        '''

        _mol = Chem.Mol(self.mol)
        _mol.UpdatePropertyCache(strict=False)
        _mol = remove_atom_map(_mol)
        args = self.args
        ps = Chem.SmilesParserParams()
        ps.sanitize = False

        if args.explicit_H:
            _mol = Chem.AddHs(_mol, explicitOnly=False)
            ps.removeHs = False

        start_idx = idx
        for atom in _mol.GetAtoms():
            atom.SetIsotope(idx)
            idx+=1

        self.atom_mapping = [i for i in range(start_idx, idx)]
        self.smiles_w_isotope = Chem.MolToSmiles(_mol)
        self.smiles_w_mapping = Chem.MolToSmiles(isotope_to_atommap(_mol))
        self.mol = Chem.MolFromSmiles(self.smiles_w_isotope, ps)
        self.identity = 'reactant'

        if args.explicit_H:
            plain = Chem.MolToSmiles(remove_atom_map(_mol, isotope=True))
            self.smiles = plain

        return idx

    def add_product(self):
        args = self.args
        _mol = Chem.Mol(self.mol)
        _mol.UpdatePropertyCache(strict=False)
        _mol = remove_atom_map(_mol)

        if args.explicit_H:
            _mol = Chem.AddHs(_mol, explicitOnly=False)
            self.mol = _mol

        self.identity = 'product'
        if args.explicit_H:
            plain = Chem.MolToSmiles(remove_atom_map(_mol, isotope=True))
            self.smiles = plain

    def add_intermediate(self):
        '''
        Product mol is isotope-labeld
        '''
        _mol = self.mol
        _mol.UpdatePropertyCache(strict=False)
        self.smiles_w_isotope = Chem.MolToSmiles(_mol)
        atom_map_mol, atom_mapping = isotope_to_atommap(_mol, get_map_list=True)
        self.smiles_w_mapping = Chem.MolToSmiles(atom_map_mol)
        self.atom_mapping = atom_mapping
        if not self.identity:
            self.identity = 'intermediate'  # Before checking if it is a recorded product, it will be intermediate.




if __name__ == '__main__':
    import argparse
    def parse_arguments():
        parser = argparse.ArgumentParser("Set the arguments for mechanistic dataset generation")
        parser.add_argument("--explicit_H", help="SMILES with Explicit H",
                            type=bool, default=True)
        return parser.parse_args()

    args = parse_arguments()

    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    #Testing the reactants
    O_smi = 'BrC1=C(Br)C=CC=C1C'
    mol = Chem.MolFromSmiles(O_smi, ps)
    O_mol = Molecule_Node(mol, args)
    i = O_mol.add_reactant(1)
    print(i)

    #Testing the intermediates
    prod_smi = '[Cl:14][Pd:15][P+:40]([c:39]1[cH:38][cH:37][cH:36][cH:54][cH:53]1)([c:41]1[cH:42][cH:43][cH:44][cH:45][cH:46]1)[c:47]1[cH:48][cH:49][cH:50][cH:51][cH:52]1'
    mol = Chem.MolFromSmiles(prod_smi,sanitize=False)
    mol = atommap_to_isotope(mol)
    p_mol = Molecule_Node(mol, args)


    p_mol.add_intermediate()

