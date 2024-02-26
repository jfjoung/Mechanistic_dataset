import os
import sys
import networkx as nx
from utils.parsing import parse_arguments
from utils.apply_mechanistic_template import get_mechanistic_network, elementary_reaction
import json

def generate_mechanism_for_one_reaction(rxn, args):
    # os.makedirs(args.save, exist_ok=True)
    print(rxn)
    G_list = get_mechanistic_network(rxn, v=False, simple=args.simple)
    elem_dict=dict()
    for G in G_list:
        try:
            elem_rxns = elementary_reaction(G, v=False, byproduct=args.byproduct, spectator=args.spectator, full=args.full,
                                            end=args.end, plain=args.plain,reagent=args.reagent)
        except: pass

    return

if __name__ == '__main__':
    args = parse_arguments()
    print(args)

    with open(args.data, 'r') as file:
        for line in file:
            rxn = json.loads(line.strip())
            generate_mechanism_for_one_reaction(rxn, args)
