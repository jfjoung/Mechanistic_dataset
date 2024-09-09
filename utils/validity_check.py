from rdkit import Chem
import numpy as np
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


bt_to_electron = {Chem.rdchem.BondType.SINGLE: 2,
                  Chem.rdchem.BondType.DOUBLE: 4,
                  Chem.rdchem.BondType.TRIPLE: 6,
                  Chem.rdchem.BondType.AROMATIC: 3}

tbl = Chem.GetPeriodicTable()


def bond_features(bond):
    bt = bond.GetBondType()

    return bt_to_electron[bt]
    # np.array([bt == Chem.rdchem.BondType.SINGLE, bt == Chem.rdchem.BondType.DOUBLE, bt == Chem.rdchem.BondType.TRIPLE, bt == Chem.rdchem.BondType.AROMATIC, bond.GetIsConjugated(), bond.IsInRing()], dtype=np.float32)


def count_lone_pairs(a):
    v = tbl.GetNOuterElecs(a.GetAtomicNum())
    c = a.GetFormalCharge()
    b = sum([bond.GetBondTypeAsDouble() for bond in a.GetBonds()])
    h = a.GetTotalNumHs()
    return v - c - b - h


def get_BE_matrix(r, kekule=False):
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
    # if not kekule:
    ps.sanitize = False
    rmol = Chem.MolFromSmiles(r, ps)
    # Chem.Kekulize(rmol)
    max_natoms = len(rmol.GetAtoms())
    f = np.zeros((max_natoms, max_natoms))

    map_dict = {}

    for atom in rmol.GetAtoms():
        lone_pair = count_lone_pairs(atom)
        f[atom.GetIntProp('molAtomMapNumber') - 1, atom.GetIntProp('molAtomMapNumber') - 1] = lone_pair
        if atom.GetIntProp('molAtomMapNumber') not in map_dict:
            map_dict[atom.GetIntProp('molAtomMapNumber')] = atom.GetAtomicNum()
        else:
            return np.zeros((max_natoms, max_natoms)), {}, max_natoms

    for bond in rmol.GetBonds():
        a1 = bond.GetBeginAtom().GetIntProp('molAtomMapNumber') - 1
        a2 = bond.GetEndAtom().GetIntProp('molAtomMapNumber') - 1
        f[(a1, a2)] = f[(a2, a1)] = bond_features(bond) / 2  # so that bond electron diff matrix sums up to 0

    return f, map_dict, max_natoms

def check_reaction_validity(smi):
    rsmi, _ , psmi = smi.split('>')

    rmat, rdict, max_natoms = get_BE_matrix(rsmi)
    pmat, pdict, _ = get_BE_matrix(psmi)


    if max_natoms != len(rmat):
        return False
    if not rdict or not pdict:
        return False
    elif rdict != pdict:
        return False
    elif sum(sum(rmat)) != sum(sum(pmat)):
        return False
    else:
        return True

def remapping(rxn_smi):
    ps = Chem.SmilesParserParams()
    ps.removeHs = False
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


if __name__ == '__main__':

    #Overall reaction
    rxn1 = '[CH3:10][B:11]([OH:12])[OH:13].[CH3:1][c:2]1[cH:3][cH:4][cH:5][c:6]([Br:7])[c:8]1[Br:9].[Cl:14][Pd:15][Cl:16].[OH2:55].[OH2:56].[OH2:57].[cH:17]1[cH:18][cH:19][c:20]([P:21]([c:22]2[cH:23][cH:24][cH:25][cH:26][cH:27]2)[c:28]2[cH:29][cH:30][cH:31][cH:32][cH:33]2)[cH:34][cH:35]1.[OH2:58].[cH:36]1[cH:37][cH:38][c:39]([P:40]([c:41]2[cH:42][cH:43][cH:44][cH:45][cH:46]2)[c:47]2[cH:48][cH:49][cH:50][cH:51][cH:52]2)[cH:53][cH:54]1>>[CH3:1][c:2]1[cH:3][cH:4][cH:5][c:6]([CH3:10])[c:8]1[Br:9].[Pd:15].[B:11]([OH:12])([OH:13])[OH:57].[BrH:7].[Cl-:16].[ClH:14].[OH3+:56].[cH:17]1[cH:18][cH:19][c:20]([P:21]([c:22]2[cH:23][cH:24][cH:25][cH:26][cH:27]2)([c:28]2[cH:29][cH:30][cH:31][cH:32][cH:33]2)=[O:55])[cH:34][cH:35]1.[OH2:58].[cH:36]1[cH:37][cH:38][c:39]([P:40]([c:41]2[cH:42][cH:43][cH:44][cH:45][cH:46]2)[c:47]2[cH:48][cH:49][cH:50][cH:51][cH:52]2)[cH:53][cH:54]1'
    print(check_reaction_validity(rxn1))

    # Elementary reaction
    rxn2 = '[C:121]1([H:122])([H:123])[C:124]([H:125])([H:126])[O:127][C:128]([H:129])([H:130])[C:131]1([H:132])[H:133].[C:1]([C:2]([C:3]([Al:4]([C:5]([C:6]([C:7]([H:8])([H:9])[H:10])([C:11]([H:12])([H:13])[H:14])[H:15])([H:16])[H:17])[H:18])([H:19])[H:20])([C:21]([H:22])([H:23])[H:24])[H:25])([H:26])([H:27])[H:28].[C:29]([Si:30]([N-:31][Si:32]([C:33]([H:34])([H:35])[H:36])([C:37]([H:38])([H:39])[H:40])[C:41]([H:42])([H:43])[H:44])([C:45]([H:46])([H:47])[H:48])[C:49]([H:50])([H:51])[H:52])([H:53])([H:54])[H:55].[Li+:56].[O:57]([C:58]([C@@:59]([C:60]([H:61])([H:62])[H:63])([C@:64]([N:65]([C:66]1([c:67]2[c:68]([H:69])[c:70]([H:71])[c:72]([H:73])[c:74]([H:75])[c:76]2[H:77])[c:78]2[c:79]([H:80])[c:81]([H:82])[c:83]([H:84])[c:85]([H:86])[c:87]2-[c:88]2[c:89]1[c:90]([H:91])[c:92]([H:93])[c:94]([H:95])[c:96]2[H:97])[H:98])([C:99](=[O:100])[O:101][C:102]([c:103]1[c:104]([H:105])[c:106]([H:107])[c:108]([H:109])[c:110]([H:111])[c:112]1[H:113])([H:114])[H:115])[H:116])[H:117])([H:118])[H:119])[H:120]>>[C:121]1([H:122])([H:123])[C:124]([H:125])([H:126])[O:127][C:128]([H:129])([H:130])[C:131]1([H:132])[H:133].[C:1]([C:2]([C:3]([Al:4]([C:5]([C:6]([C:7]([H:8])([H:9])[H:10])([C:11]([H:12])([H:13])[H:14])[H:15])([H:16])[H:17])[H:18])([H:19])[H:20])([C:21]([H:22])([H:23])[H:24])[H:25])([H:26])([H:27])[H:28].[C:29]([Si:30]([N-:31][Si:32]([C:33]([H:34])([H:35])[H:36])([C:37]([H:38])([H:39])[H:40])[C:41]([H:42])([H:43])[H:44])([C:45]([H:46])([H:47])[H:48])[C:49]([H:50])([H:51])[H:52])([H:53])([H:54])[H:55].[Li+:56].[O+:57]1([H:120])[C:58]([H:118])([H:119])[C@:59]([C:60]([H:61])([H:62])[H:63])([H:117])[C@@:64]([N:65]([C:66]2([c:67]3[c:68]([H:69])[c:70]([H:71])[c:72]([H:73])[c:74]([H:75])[c:76]3[H:77])[c:78]3[c:79]([H:80])[c:81]([H:82])[c:83]([H:84])[c:85]([H:86])[c:87]3-[c:88]3[c:89]2[c:90]([H:91])[c:92]([H:93])[c:94]([H:95])[c:96]3[H:97])[H:98])([H:116])[C:99]1([O-:100])[O:101][C:102]([c:103]1[c:104]([H:105])[c:106]([H:107])[c:108]([H:109])[c:110]([H:111])[c:112]1[H:113])([H:114])[H:115]'
    print(check_reaction_validity(rxn2))

    # Elementary reaction, atom map collision
    rxn3 = '[C:1]1([H:2])([H:3])[C:4]([H:5])([H:6])[O:7][C:8]([H:9])([H:10])[C:11]1([H:12])[H:13].[C:1]([C:2]([C:3]([Al:4]([C:5]([C:6]([C:7]([H:8])([H:9])[H:10])([C:11]([H:12])([H:13])[H:14])[H:15])([H:16])[H:17])[H:18])([H:19])[H:20])([C:21]([H:22])([H:23])[H:24])[H:25])([H:26])([H:27])[H:28].[C:29]([Si:30]([N-:31][Si:32]([C:33]([H:34])([H:35])[H:36])([C:37]([H:38])([H:39])[H:40])[C:41]([H:42])([H:43])[H:44])([C:45]([H:46])([H:47])[H:48])[C:49]([H:50])([H:51])[H:52])([H:53])([H:54])[H:55].[Li+:56].[O:57]([C:58]([C@@:59]([C:60]([H:61])([H:62])[H:63])([C@:64]([N:65]([C:66]1([c:67]2[c:68]([H:69])[c:70]([H:71])[c:72]([H:73])[c:74]([H:75])[c:76]2[H:77])[c:78]2[c:79]([H:80])[c:81]([H:82])[c:83]([H:84])[c:85]([H:86])[c:87]2-[c:88]2[c:89]1[c:90]([H:91])[c:92]([H:93])[c:94]([H:95])[c:96]2[H:97])[H:98])([C:99](=[O:100])[O:101][C:102]([c:103]1[c:104]([H:105])[c:106]([H:107])[c:108]([H:109])[c:110]([H:111])[c:112]1[H:113])([H:114])[H:115])[H:116])[H:117])([H:118])[H:119])[H:120]>>[C:1]1([H:2])([H:3])[C:4]([H:5])([H:6])[O:7][C:8]([H:9])([H:10])[C:11]1([H:12])[H:13].[C:1]([C:2]([C:3]([Al:4]([C:5]([C:6]([C:7]([H:8])([H:9])[H:10])([C:11]([H:12])([H:13])[H:14])[H:15])([H:16])[H:17])[H:18])([H:19])[H:20])([C:21]([H:22])([H:23])[H:24])[H:25])([H:26])([H:27])[H:28].[C:29]([Si:30]([N-:31][Si:32]([C:33]([H:34])([H:35])[H:36])([C:37]([H:38])([H:39])[H:40])[C:41]([H:42])([H:43])[H:44])([C:45]([H:46])([H:47])[H:48])[C:49]([H:50])([H:51])[H:52])([H:53])([H:54])[H:55].[Li+:56].[O+:57]1([H:120])[C:58]([H:118])([H:119])[C@:59]([C:60]([H:61])([H:62])[H:63])([H:117])[C@@:64]([N:65]([C:66]2([c:67]3[c:68]([H:69])[c:70]([H:71])[c:72]([H:73])[c:74]([H:75])[c:76]3[H:77])[c:78]3[c:79]([H:80])[c:81]([H:82])[c:83]([H:84])[c:85]([H:86])[c:87]3-[c:88]3[c:89]2[c:90]([H:91])[c:92]([H:93])[c:94]([H:95])[c:96]3[H:97])[H:98])([H:116])[C:99]1([O-:100])[O:101][C:102]([c:103]1[c:104]([H:105])[c:106]([H:107])[c:108]([H:109])[c:110]([H:111])[c:112]1[H:113])([H:114])[H:115]'
    print(check_reaction_validity(rxn3))

    # Elementary reaction, missing atom composition
    rxn4 = '[C:1]([C:2]([C:3]([Al:4]([C:5]([C:6]([C:7]([H:8])([H:9])[H:10])([C:11]([H:12])([H:13])[H:14])[H:15])([H:16])[H:17])[H:18])([H:19])[H:20])([C:21]([H:22])([H:23])[H:24])[H:25])([H:26])([H:27])[H:28].[C:29]([Si:30]([N-:31][Si:32]([C:33]([H:34])([H:35])[H:36])([C:37]([H:38])([H:39])[H:40])[C:41]([H:42])([H:43])[H:44])([C:45]([H:46])([H:47])[H:48])[C:49]([H:50])([H:51])[H:52])([H:53])([H:54])[H:55].[Li+:56].[O:57]([C:58]([C@@:59]([C:60]([H:61])([H:62])[H:63])([C@:64]([N:65]([C:66]1([c:67]2[c:68]([H:69])[c:70]([H:71])[c:72]([H:73])[c:74]([H:75])[c:76]2[H:77])[c:78]2[c:79]([H:80])[c:81]([H:82])[c:83]([H:84])[c:85]([H:86])[c:87]2-[c:88]2[c:89]1[c:90]([H:91])[c:92]([H:93])[c:94]([H:95])[c:96]2[H:97])[H:98])([C:99](=[O:100])[O:101][C:102]([c:103]1[c:104]([H:105])[c:106]([H:107])[c:108]([H:109])[c:110]([H:111])[c:112]1[H:113])([H:114])[H:115])[H:116])[H:117])([H:118])[H:119])[H:120]>>[C:121]1([H:122])([H:123])[C:124]([H:125])([H:126])[O:127][C:128]([H:129])([H:130])[C:131]1([H:132])[H:133].[C:1]([C:2]([C:3]([Al:4]([C:5]([C:6]([C:7]([H:8])([H:9])[H:10])([C:11]([H:12])([H:13])[H:14])[H:15])([H:16])[H:17])[H:18])([H:19])[H:20])([C:21]([H:22])([H:23])[H:24])[H:25])([H:26])([H:27])[H:28].[C:29]([Si:30]([N-:31][Si:32]([C:33]([H:34])([H:35])[H:36])([C:37]([H:38])([H:39])[H:40])[C:41]([H:42])([H:43])[H:44])([C:45]([H:46])([H:47])[H:48])[C:49]([H:50])([H:51])[H:52])([H:53])([H:54])[H:55].[Li+:56].[O+:57]1([H:120])[C:58]([H:118])([H:119])[C@:59]([C:60]([H:61])([H:62])[H:63])([H:117])[C@@:64]([N:65]([C:66]2([c:67]3[c:68]([H:69])[c:70]([H:71])[c:72]([H:73])[c:74]([H:75])[c:76]3[H:77])[c:78]3[c:79]([H:80])[c:81]([H:82])[c:83]([H:84])[c:85]([H:86])[c:87]3-[c:88]3[c:89]2[c:90]([H:91])[c:92]([H:93])[c:94]([H:95])[c:96]3[H:97])[H:98])([H:116])[C:99]1([O-:100])[O:101][C:102]([c:103]1[c:104]([H:105])[c:106]([H:107])[c:108]([H:109])[c:110]([H:111])[c:112]1[H:113])([H:114])[H:115]'
    print(check_reaction_validity(rxn4))

    # Elementary reaction, different atom composision
    rxn5 = '[Si:121]1([F:122])([Cl:123])[C:124]([H:125])([H:126])[O:127][C:128]([H:129])([H:130])[C:131]1([H:132])[H:133].[C:1]([C:2]([C:3]([Al:4]([C:5]([C:6]([C:7]([H:8])([H:9])[H:10])([C:11]([H:12])([H:13])[H:14])[H:15])([H:16])[H:17])[H:18])([H:19])[H:20])([C:21]([H:22])([H:23])[H:24])[H:25])([H:26])([H:27])[H:28].[C:29]([Si:30]([N-:31][Si:32]([C:33]([H:34])([H:35])[H:36])([C:37]([H:38])([H:39])[H:40])[C:41]([H:42])([H:43])[H:44])([C:45]([H:46])([H:47])[H:48])[C:49]([H:50])([H:51])[H:52])([H:53])([H:54])[H:55].[Li+:56].[O:57]([C:58]([C@@:59]([C:60]([H:61])([H:62])[H:63])([C@:64]([N:65]([C:66]1([c:67]2[c:68]([H:69])[c:70]([H:71])[c:72]([H:73])[c:74]([H:75])[c:76]2[H:77])[c:78]2[c:79]([H:80])[c:81]([H:82])[c:83]([H:84])[c:85]([H:86])[c:87]2-[c:88]2[c:89]1[c:90]([H:91])[c:92]([H:93])[c:94]([H:95])[c:96]2[H:97])[H:98])([C:99](=[O:100])[O:101][C:102]([c:103]1[c:104]([H:105])[c:106]([H:107])[c:108]([H:109])[c:110]([H:111])[c:112]1[H:113])([H:114])[H:115])[H:116])[H:117])([H:118])[H:119])[H:120]>>[C:121]1([H:122])([H:123])[C:124]([H:125])([H:126])[O:127][C:128]([H:129])([H:130])[C:131]1([H:132])[H:133].[C:1]([C:2]([C:3]([Al:4]([C:5]([C:6]([C:7]([H:8])([H:9])[H:10])([C:11]([H:12])([H:13])[H:14])[H:15])([H:16])[H:17])[H:18])([H:19])[H:20])([C:21]([H:22])([H:23])[H:24])[H:25])([H:26])([H:27])[H:28].[C:29]([Si:30]([N-:31][Si:32]([C:33]([H:34])([H:35])[H:36])([C:37]([H:38])([H:39])[H:40])[C:41]([H:42])([H:43])[H:44])([C:45]([H:46])([H:47])[H:48])[C:49]([H:50])([H:51])[H:52])([H:53])([H:54])[H:55].[Li+:56].[O+:57]1([H:120])[C:58]([H:118])([H:119])[C@:59]([C:60]([H:61])([H:62])[H:63])([H:117])[C@@:64]([N:65]([C:66]2([c:67]3[c:68]([H:69])[c:70]([H:71])[c:72]([H:73])[c:74]([H:75])[c:76]3[H:77])[c:78]3[c:79]([H:80])[c:81]([H:82])[c:83]([H:84])[c:85]([H:86])[c:87]3-[c:88]3[c:89]2[c:90]([H:91])[c:92]([H:93])[c:94]([H:95])[c:96]3[H:97])[H:98])([H:116])[C:99]1([O-:100])[O:101][C:102]([c:103]1[c:104]([H:105])[c:106]([H:107])[c:108]([H:109])[c:110]([H:111])[c:112]1[H:113])([H:114])[H:115]'
    print(check_reaction_validity(rxn5))