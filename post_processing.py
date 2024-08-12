from utils.explicit_H import modify_explicit_H
from utils.validity_check import check_reaction_validity, remapping
from scripts.setting import Args
import os
from tqdm import tqdm


def main(args):

    with open(args.save, 'r') as file:
        lines = file.readlines()

    base_file_root, _ = os.path.splitext(args.save)
    postprocess_file_path = f"{base_file_root}_explicit.txt"

    num_rxn = 0

    with open(postprocess_file_path, 'w') as fout:
        for implicit_rxn in tqdm(lines):
            H_rxn, _ = modify_explicit_H(implicit_rxn)
            if check_reaction_validity(H_rxn):
                fout.write(f'{remapping(H_rxn)}\n')
                num_rxn += 1
                

if __name__ == '__main__':
    args=Args()
    main(args)