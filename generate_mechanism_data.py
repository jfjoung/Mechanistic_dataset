import argparse
from scripts.generate_mech_data import generate_mechdata_known_condition, generate_mechdata_unknown_condition

def parse_arguments():
    parser = argparse.ArgumentParser("Set the arguments for mechanistic dataset generation")

    parser.add_argument("--all_info", help="Save all the information you can get, or you can get only elementary steps if not",
                        type=bool, default=False)

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
    parser.add_argument("--simple", help="Use all simple path finding instead of shortest path during reaction network generation",
                        type=bool, default=False)


    parser.add_argument('--data', help='Path to the reaction data',
                        type=str, default='./data/test_data.txt')
    parser.add_argument('--save', help='Path to the saving data',
                        type=str, default='./results/elementary_reaction.txt')


    #TODO: set the verbosity level
    parser.add_argument("--verbosity", help="control the verbosity",
                       type=int, default=1)

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    print(args)
    '''
    If you have a list of reactions in the form of 
    {'reaction_name': NameRXN name,
    'reaction_smiles': Reaction SMILES,
    'conditions': a list of conditions for a given reaction name},
    then RUN generate_mechdata_known_condition(args)
    
    or
    if you have a reaction string of 'reaction_smiles NameRXN_name', 
    then RUN generate_mechdata_unknown_condition(args)
    '''

    # generate_mechdata_known_condition(args)
    generate_mechdata_unknown_condition(args)
