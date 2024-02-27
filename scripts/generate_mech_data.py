import os
import sys
sys.path.append("../")
sys.path.append("../utils")
import networkx as nx
from tqdm import tqdm
from utils.parsing import parse_arguments
from utils.apply_mechanistic_template import get_mechanistic_network, elementary_reaction, flatten_list
import json
# import signal

def generate_mechanism_for_one_reaction(rxn, args):
    G_dict = get_mechanistic_network(rxn, v=False, simple=args.simple)
    if args.all_info:
        elem_dict = dict()
    else:
        elem_list = list()

    for cond, G in G_dict.items():
        try:
            elem_rxns = elementary_reaction(G, v=False, byproduct=args.byproduct, spectator=args.spectator, full=args.full,
                                            end=args.end, plain=args.plain,reagent=args.reagent)
            if args.all_info:
                elem_dict[cond]={'Reaction graph': nx.node_link_data(G),
                                 'Elementary steps': elem_rxns}
            else:
                elem_list.append(elem_rxns)

        except Exception as e: pass

    if args.all_info:
        rxn['Mechanism'] = elem_dict
        return rxn
    else:
        return flatten_list(elem_list)

def generate_batch(args):
    with open(args.data, 'r') as file, open(args.save, 'w') as fout:
        for i, line in tqdm(enumerate(file)):
            # if i<3604:
            #     continue
            # print(i)
            # print(line)
            try:
                rxn = json.loads(line.strip())
                new_rxn = generate_mechanism_for_one_reaction(rxn, args)
                if args.all_info:
                    fout.write('{}\n'.format(new_rxn))
                else:
                    for step_rxn in new_rxn:
                        fout.write('{}\n'.format(step_rxn))
            except: pass
