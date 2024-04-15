import argparse
from scripts.analysis import analysis, mol_analysis, ERS_analysis, template_analysis


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False

def parse_arguments():
    parser = argparse.ArgumentParser("Set the arguments for analysis of mechanistic dataset")
    parser.add_argument('--data', help='Path to the reaction data to be analyzed, it should be jsonl format',
                        type=str, default='./results/uspto_all.jsonl')
    parser.add_argument('--weight_bin', help='Bins for molecular weight',
                        type=int, default=20)
    parser.add_argument('--atom_bin', help='Bins for heavy atom',
                        type=int, default=2)
    parser.add_argument("--process", help="The number of worker processes",
                        type=int, default=30)
    parser.add_argument("--log_scale", help="Draw figure in log scale",
                        type=str2bool, default=True)
    return parser.parse_args()

if __name__ == '__main__':
    args =  parse_arguments()
    # analysis(args)
    # mol_analysis(args)
    ERS_analysis(args)
    # template_analysis(args)