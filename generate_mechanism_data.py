import os
import argparse
import logging
from datetime import datetime
from scripts.generate_mech_data import generate_mechdata

def parse_arguments():
    parser = argparse.ArgumentParser("Set the arguments for mechanistic dataset generation")

    # Arguments for applying reaction templates to get every elementary reactions
    parser.add_argument("--num_combination", help="The number of possible reactions to be considered, it prevents combinatorial explosions",
                        type=int, default=12)
    parser.add_argument("--uni_rxn", help="Allow the unimolecular reactions",
                        type=bool, default=True)
    parser.add_argument("--proton", help="Allow the proton-balanced reactions",
                        type=bool, default=True)
    parser.add_argument("--stoichiometry", help="Duplicate the reactants when needed",
                        type=bool, default=True)

    # Arguments for handling the reaction network
    parser.add_argument("--simple", help="Use all simple path finding instead of shortest path during reaction network generation",
                        type=bool, default=False)
    parser.add_argument("--do_not_pruning", help="Get full reaction network instead of the pruned network containing only paths toward product",
                        type=bool, default=False)

    # Arguments for extracting reaction SMARTS from the reaction network
    parser.add_argument("--byproduct", help="Add produced byproduct in reaction SMARTS",
                       type=bool, default=False)
    parser.add_argument("--spectator", help="Add spectator in reaction SMARTS",
                       type=bool, default=False)
    parser.add_argument("--full", help="Get overall reaction instead of elementary steps",
                       type=bool, default=False)
    parser.add_argument("--end", help="Get termination reactions (products>>products)",
                       type=bool, default=False)
    parser.add_argument("--plain", help="Get reaction SMARTS without atom-mapping",
                       type=bool, default=False)
    parser.add_argument("--reagent", help="Locate reagents at middle of the reaction (reactants>reagents>products) instead of both sides of reactants and products.",
                       type=bool, default=False)

    parser.add_argument("--num_cycles", help="The maximum number of cycles allowed in the reaction graph",
                       type=int, default=9)
    parser.add_argument("--num_reaction_node", help="The maximum number of reaction nodes allowed in the reaction graph",
                       type=int, default=50)

    # Arguments for data loading and saving
    parser.add_argument('--data', help='Path to the reaction data',
                        type=str, default='./data/test_data.txt')
    parser.add_argument('--save', help='Path to the saving data',
                        type=str, default='./results/elementary_reaction.txt')
    parser.add_argument("--all_info", help="Save all the information you can get, or you can get only elementary steps if not",
                        type=bool, default=False)
    parser.add_argument("--stat", help="Save the statistics dictionary",
                        type=bool, default=True)

    parser.add_argument("--verbosity", help="control the verbosity; 0: silent, 1: prints the critical errors, 2: prints some details, 3: prints all",
                       type=int, default=0)

    return parser.parse_args()

if __name__ == '__main__':
    '''
    You need a reaction string of 'reaction_smiles NameRXN_name', 
    generate_mechdata_unknown_condition(args)
    '''
    args = parse_arguments()
    os.makedirs('./logs/', exist_ok=True)
    time = datetime.now().strftime('%Y-%m-%d_%H-%M')
    logging.basicConfig(filename=f'./logs/generate_mechanism_data_{time}.log', level=logging.DEBUG, format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.info('Starting mechanistic dataset generation...')
    logging.info(f'Arguments')
    for key, value in vars(args).items():
        logging.info('{}: {}'.format(key, value))
    generate_mechdata(args)
