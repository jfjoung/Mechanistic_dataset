import os
import argparse
import logging
from datetime import datetime
from scripts.generate_mech_data import generate_mechdata_multiprocess

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False

def parse_arguments():
    parser = argparse.ArgumentParser("Set the arguments for mechanistic dataset generation")

    # Arguments for applying reaction templates to get every elementary reactions
    parser.add_argument("--num_combination", help="The number of possible reactions to be considered, it prevents combinatorial explosions",
                        type=int, default=12)
    parser.add_argument("--uni_rxn", help="Allow the unimolecular reactions",
                        type=str2bool, default=True)
    parser.add_argument("--proton", help="Allow the proton-balanced reactions",
                        type=str2bool, default=True)
    parser.add_argument("--max_num_temp", help="The maximum number of templates allowed to consider. If uni_rxn or proton is true, the combinatorial explosion could occur",
                        type=int, default=50)
    parser.add_argument("--stoichiometry", help="Duplicate the reactants when needed",
                        type=str2bool, default=True)

    # Arguments for handling the reaction network
    parser.add_argument("--simple", help="Use all simple path finding instead of shortest path during reaction network generation",
                        type=str2bool, default=False)
    parser.add_argument("--do_not_pruning", help="Get full reaction network instead of the pruned network containing only paths toward product",
                        type=str2bool, default=False)

    # Arguments for preventing the combinatorial explosions in the reaction network
    parser.add_argument("--num_cycles", help="The maximum number of cycles allowed in the reaction graph",
                       type=int, default=9)
    parser.add_argument("--num_reaction_node", help="The maximum number of reaction nodes allowed in the reaction graph",
                       type=int, default=50)

    # Arguments for extracting reaction SMARTS from the reaction network
    '''
    spectator: It will add spectators for all reactions. Reactants that have not yet participated in the reaction are also included. 
               e.g. Reactants.Spectators>>Products.Spectators
    byproduct: It will add produced byproducts as spectators for downstream reaction. 
               e.g. upstream reaction; reactants>>intermediates.byproducts
                    downstream reaction; reactants.intermediates.byproducts>>new_intermediates.byproducts
    full: It returns the overall reactions.
    end: It returns the termination reactions. e.g. products>>products
    reagent: It changes the SMARTS format from reactants.spectators>>products.spectators to reactants>spectators>products
    '''

    parser.add_argument("--byproduct", help="Add produced byproduct in reaction SMARTS",
                       type=str2bool, default=False)
    parser.add_argument("--spectator", help="Add spectator in reaction SMARTS",
                       type=str2bool, default=False)
    parser.add_argument("--full", help="Get overall reaction instead of elementary steps",
                       type=str2bool, default=False)
    parser.add_argument("--end", help="Get termination reactions (products>>products)",
                       type=str2bool, default=False)
    parser.add_argument("--plain", help="Get reaction SMARTS without atom-mapping",
                       type=str2bool, default=False)
    parser.add_argument("--reagent", help="Locate reagents at middle of the reaction (reactants>reagents>products) instead of both sides of reactants and products.",
                       type=str2bool, default=False)
    parser.add_argument("--remapping", help="Re-atom-mapping for each elementary reaction to start with 1.",
                       type=str2bool, default=False)
    # Arguments for data loading and saving
    parser.add_argument('--data', help='Path to the reaction data',
                        type=str, default='./data/test_data.txt')
    parser.add_argument('--save', help='Path to the saving data',
                        type=str, default='./results/test.txt')
    parser.add_argument("--all_info", help="Save all the information you can get, or you can get only elementary steps if not",
                        type=str2bool, default=False)
    parser.add_argument("--stat", help="Save the statistics dictionary",
                        type=str2bool, default=True)
    parser.add_argument("--rxn_class", help="Specify the reaction class to work on. If not specified, it will work on all classes.",
                        type=str, default='')
    parser.add_argument("--process", help="The number of worker processes",
                        type=int, default=10)
    parser.add_argument("--verbosity", help="control the verbosity; 0: silent, 1: prints the critical errors, 2: prints some details (do not recommend more than 10 reactions), 3: prints the process, 4: prints all",
                       type=int, default=0)

    return parser.parse_args()

if __name__ == '__main__':
    '''
    You need a reaction string of 'reaction_smiles NameRXN_name', 
    generate_mechdata_multiprocess(args)
    '''
    args = parse_arguments()
    os.makedirs('./logs/', exist_ok=True)
    time = datetime.now().strftime('%Y-%m-%d_%H-%M')
    logging.basicConfig(filename=f'./logs/generate_mechanism_data_{time}.log', level=logging.DEBUG, format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.info('Starting mechanistic dataset generation...')
    logging.info(f'Arguments')
    for key, value in vars(args).items():
        logging.info('{}: {}'.format(key, value))
    generate_mechdata_multiprocess(args)