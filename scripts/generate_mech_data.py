import sys
sys.path.append('/mnt/home/jjoung/Mechanistic_dataset_generation')

from utils.parsing import parse_arguments
from utils.apply_mechanistic_template import get_mechanistic_network, elementary_reaction
import json

def generate_mechanism_for_one_reaction(rxn, args):
    G_list = get_mechanistic_network(rxn, v=False, simple=args.simple)
    for G in G_list:
        elem_list = elementary_reaction(G, v=False, byproduct=args.byproduct, spectator=args.spectator, full=args.full,
                                        end=args.end, plain=args.plain,reagent=args.reagent)
        for elem_rxn in elem_list:
            print(elem_rxn)


if __name__ == '__main__':
    args = parse_arguments()
    print(args)

    with open(args.data, 'r') as file:
        for line in file:
            json_object = json.loads(line.strip())
            print(json_object)
            break