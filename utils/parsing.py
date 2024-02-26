import argparse


def parse_arguments():
    parser = argparse.ArgumentParser("Set the arguments for mechanistic dataset generation")
    parser.add_argument("--byproduct", help="Add byproduct in reaction SMARTS",
                       type=bool, default=True)
    parser.add_argument("--spectator", help="Add spectator in reaction SMARTS",
                       type=bool, default=True)
    parser.add_argument("--full", help="Get overall reaction instead of elementary steps",
                       type=bool, default=False)
    parser.add_argument("--end", help="Get termination reactions (products>>products)",
                       type=bool, default=True)
    parser.add_argument("--plain", help="Get reaction SMARTS without atom-mapping",
                       type=bool, default=False)
    parser.add_argument("--reagent", help="Locate reagents at middle of the reaction (reactants>reagents>products) instead of both sides of reactants and products.",
                       type=bool, default=False)

    parser.add_argument("--simple", help="Use all simple path finding instead of shortest path during reaction network generation",
                        type=bool, default=False)


    parser.add_argument('--data', help='Path to the reaction data',
                        type=str, default='./data/reaction_test.jsonl')
    parser.add_argument('--save', help='Path to the saving data',
                        type=str, default='./data/elementary_reaction.jsonl')


    #TODO: set the verbosity level
    parser.add_argument("--verbosity", help="control the verbosity",
                       type=int, default=0)

    return parser.parse_args()