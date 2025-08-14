# -*- coding: utf-8 -*-
"""
Deduplicate reactions by SMILES and save as a single Parquet file.
Each unique reaction gets a running Reaction_ID.
Patent_IDs are stored as lists (deduplicated) for fast querying.
"""
import os
import argparse
import pandas as pd
import logging
import psutil
import pyarrow as pa
import pyarrow.parquet as pq
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdchem
RDLogger.DisableLog('rdApp.*')

from pathlib import Path

def setup_logger(name, log_dir):
    """Set up logging configuration."""
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"{name}.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    logging.info(f"Logging to {log_file}")

metals = [
            3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
            31, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 55, 56,
            57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73,
            74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 87, 88, 89, 90, 91, 92, 93,
            94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108,
            113, 114, 115, 116
]
max_allowed_valences = {
    7 : 3, 
    8 : 2,
}


def fullSanitize (smi):
    '''
    Takes molecule smi -> fixes valence errors associated with charge/coordination by using dative bonds
    '''
    if 'Cu' in smi or 'Pd' in smi or 'Ni' in smi:
        try:
            # print(f"Sanitizing SMILES: {smi}")
            settings = Chem.SmilesParserParams()
            settings.removeHs = False
            settings.sanitize = False
            mol = Chem.RWMol(Chem.MolFromSmiles(smi, settings))
            for atom in mol.GetAtoms():
                Z = atom.GetAtomicNum()
                valence = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
                max_val = max_allowed_valences.get(Z, None)
                if max_val is not None and valence > max_val or atom.GetFormalCharge() < 0:
                    for bond in atom.GetBonds():
                        neighbor = bond.GetOtherAtom(atom)
                        if neighbor.GetAtomicNum() in metals:
                            begin_idx = atom.GetIdx()
                            end_idx = neighbor.GetIdx()
                            mol.RemoveBond(begin_idx, end_idx)
                            mol.AddBond(begin_idx, end_idx, rdchem.BondType.DATIVE)
                if atom.GetAtomicNum() == 29:
                    for bond in atom.GetBonds():
                        neighbor = bond.GetOtherAtom(atom)
                        if neighbor.GetAtomicNum() == 7 and neighbor.GetIsAromatic():
                            begin_idx = neighbor.GetIdx()
                            end_idx = atom.GetIdx()
                            mol.RemoveBond(begin_idx, end_idx)
                            mol.AddBond(begin_idx, end_idx, rdchem.BondType.DATIVE)
            quit = False
            while not quit:
                try:
                    Chem.SanitizeMol(mol)
                    return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
                except Exception as e:
                    try:
                        Chem.SanitizeMol(mol)
                        return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
                    except Exception:
                        dative_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DATIVE]
                        
                        if not dative_bonds:
                            # No more dative bonds to fix — give up
                            #print(e)
                            return None
                        
                        for bond in dative_bonds:
                            begin = bond.GetBeginAtomIdx()
                            end = bond.GetEndAtomIdx()

                            # Temporarily replace dative bond with single bond
                            mol.RemoveBond(begin, end)
                            mol.AddBond(begin, end, Chem.rdchem.BondType.SINGLE)

                            try:
                                Chem.SanitizeMol(mol)
                                # Success! Keep this change and return
                                return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
                            except Exception:
                                # Failed, revert to dative bond
                                mol.RemoveBond(begin, end)
                                mol.AddBond(begin, end, Chem.rdchem.BondType.DATIVE)

                        # If none of the single-bond replacements worked, give up
                        quit = True

        except Exception as e:
            return
            #print(f"Error sanitizing SMILES {Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)}: {e}")
    try: 
        ps = Chem.SmilesParserParams()
        ps.removeHs = False
        ps.sanitize = True
        smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi, ps), isomericSmiles=False, canonical=True)
        return smi
    except rdchem.AtomValenceException: 
        try:
            # print(f"Sanitizing SMILES: {smi}")
            settings = Chem.SmilesParserParams()
            settings.removeHs = False
            settings.sanitize = False
            mol = Chem.RWMol(Chem.MolFromSmiles(smi, settings))
            for atom in mol.GetAtoms():
                Z = atom.GetAtomicNum()
                valence = sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds())
                max_val = max_allowed_valences.get(Z, None)
                if max_val is not None and valence > max_val or atom.GetFormalCharge() < 0:
                    for bond in atom.GetBonds():
                        neighbor = bond.GetOtherAtom(atom)
                        if neighbor.GetAtomicNum() in metals:

                            begin_idx = atom.GetIdx()
                            end_idx = neighbor.GetIdx()
                            mol.RemoveBond(begin_idx, end_idx)
                            mol.AddBond(begin_idx, end_idx, rdchem.BondType.DATIVE)

            Chem.SanitizeMol(mol)
            return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)

        except Exception as e:
            print(f"Error sanitizing SMILES {smi}: {e}")
    except Exception as e:
        print(f"Error processing SMILES {smi}: {e}")

settings = Chem.SmilesParserParams()
settings.removeHs = False
settings.sanitize = False

def canonicalize_reaction_smiles(rxn_smiles):
    # Split reaction into reactants, agents, products
    parts = rxn_smiles.split(">")
    if len(parts) == 2:
        reactants, products = parts[0], parts[1]
        agents = ""
    elif len(parts) == 3:
        reactants, agents, products = parts
    else:
        raise ValueError("Invalid reaction SMILES format")

    # Canonicalize each side
    def canonicalize_side(side):
        mols = []
        for smi in side.split('.'):
            if not smi:
                continue
            sanitized = fullSanitize(smi)
            if sanitized is None:
                return None  # Fail if sanitization fails
            mol = Chem.MolFromSmiles(sanitized, settings)
            if mol is None:
                return None  # Fail if any molecule is invalid
            mols.append(mol)
        can_smiles = sorted([Chem.MolToSmiles(mol, canonical=True, allHsExplicit=True) for mol in mols])
        return ".".join(can_smiles)

    can_reactants = canonicalize_side(reactants)
    if can_reactants is None:
        return None
    can_agents = canonicalize_side(agents) if agents else ""
    if agents and can_agents is None:
        return None
    can_products = canonicalize_side(products)
    if can_products is None:
        return None

    # Rebuild canonical reaction SMILES
    if can_agents:
        return f"{can_reactants}>{can_agents}>{can_products}"
    else:
        return f"{can_reactants}>>{can_products}"
    
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

def process_line(line):
    line = line.strip()
    if not line:
        return None
    parts = line.split('|')
    if len(parts) < 2:
        return None
    SMILES = parts[0]
    # rxn_class = parts[1]
    # template_name = parts[2]
    # step = parts[3]
    # num = parts[4]
    num = parts[-1]
    can_smiles = canonicalize_reaction_smiles(SMILES)
    if can_smiles is None:
        return None
    #return (can_smiles, (can_smiles, rxn_class, template_name, step, num))
    return (can_smiles, (can_smiles, num))


def main(args):
    setup_logger("reaction_deduplication", args.log_dir)
    data_dir = Path(args.data_dir)
    file_path = data_dir / "templated.txt"
    # Change output_file to a .txt extension by default or allow user to specify
    output_file = Path(args.output_file)
    if output_file.suffix != '.txt':
        output_file = output_file.with_suffix('.txt') # Ensure it's a .txt file
        logging.warning(f"Output file extension changed to .txt: {output_file}")

    logging.info(f"Processing {file_path}")


    reactions = set()
    total_written_lines = 0
    with open(file_path, 'r', encoding='utf-8') as infile, open(output_file, 'w', encoding='utf-8') as outfile:
        with Pool(processes=cpu_count()) as pool:
            # tqdm without total, since we don't know it in advance
            for result in tqdm(pool.imap_unordered(process_line, infile, chunksize=100), desc="Deduplicating reactions (multiprocess, streaming)"):
                if result is None:
                    continue
                can_smiles, data = result
                if can_smiles not in reactions:
                    reactions.add(can_smiles)
                    outfile.write('|'.join(data) + '\n')
                    total_written_lines += 1
                    if total_written_lines % 100000 == 0:
                        logging.info(f"Wrote {total_written_lines} lines to {output_file}. Memory: {psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024:.2f} MB")

    logging.info(f"Memory used after parsing (MB): {psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024:.2f}")
    logging.info(f"Found {len(reactions)} unique reactions.")
    logging.info(f"Finished writing. Total {total_written_lines} deduplicated reactions saved to {output_file}")
    logging.info(f"Final memory used (MB): {psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024:.2f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Deduplicate reaction data and save as a text file.")
    parser.add_argument("--data_dir", type=str, default="./data/reactions", help="Input reaction data directory")
    # Change default output file name to reflect .txt output
    parser.add_argument("--output_file", type=str, default="./data/reactions_deduplicated.txt", help="Output text file")
    parser.add_argument("--log_dir", type=str, default="./logs/dataset_curation", help="Log directory")
    args = parser.parse_args()
    main(args)
