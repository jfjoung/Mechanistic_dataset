import os
from multiprocessing import Pool
from tqdm import tqdm
from utils.explicit_H import modify_explicit_H
from utils.validity_check import check_reaction_validity
from scripts.setting import Args

def process_reaction(rxn):
    """Process a single reaction, modify explicit H and check validity."""
    args = Args()
    # try:
    #     H_rxn, _ = modify_explicit_H(rxn)
    # except:
    #     return None

    if '|' in rxn:
        new_rxn, _ = rxn.split('|')
    else:
        new_rxn = rxn
    try:
        if check_reaction_validity(new_rxn, args.sanitize):
            return rxn # remapping(H_rxn)
        else: return None
    except:
        return None

def save_valid_reactions(file_path, valid_reactions):
    """Save valid reactions to the specified file."""
    with open(file_path, 'w') as fout:
        for rxn in tqdm(valid_reactions, total=len(valid_reactions), desc="Saving reactions"):
            if rxn:  # Ensure None values are not written
                fout.write(f'{rxn}')

def main(args):
    try:
        with open(args.save, 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        print(f"File not found: {args.save}")
        return

    base_file_root, _ = os.path.splitext(args.save)
    postprocess_file_path = f"{base_file_root}_explicit.txt"

    # Using multiprocessing to process reactions in parallel
    with Pool(args.process) as pool:
        valid_reactions = list(tqdm(pool.imap(process_reaction, lines), total=len(lines), desc="Processing reactions"))

    # Filter out None values that indicate invalid reactions
    valid_reactions = [rxn for rxn in valid_reactions if rxn is not None]

    save_valid_reactions(postprocess_file_path, valid_reactions)

    num_rxn = len(valid_reactions)
    percent_remaining = (num_rxn / len(lines)) * 100
    print(f'Total {num_rxn} reactions remain out of {len(lines)} ({percent_remaining:.1f}%)')

if __name__ == '__main__':
    args = Args()
    main(args)

