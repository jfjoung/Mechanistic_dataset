import os
import sys
import logging
import networkx as nx
from tqdm import tqdm
import json
from utils.apply_mechanistic_template import get_mechanistic_network, elementary_reaction, flatten_list, reagent_matching_for_single_reaction


def generate_mechanism_for_one_reaction(rxn, args):
    logging.info(f'Generating mechanism for {rxn['reaction_smiles']}...')
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

def generate_mechdata_known_condition(args):
    '''
    Input format is
    {'reaction_name': NameRXN name,
    'reaction_smiles': Reaction SMILES,
    'conditions': a list of conditions for a given reaction name}
    '''

    # Make directory for saving file
    save_dir_path=os.path.dirname(args.save)
    os.makedirs(save_dir_path, exist_ok=True)

    with open(args.data, 'r') as file, open(args.save, 'w') as fout:
        for i, line in tqdm(enumerate(file)):
            try:
                rxn = json.loads(line.strip())
                new_rxn = generate_mechanism_for_one_reaction(rxn, args)
                if args.all_info:
                    fout.write('{}\n'.format(new_rxn))
                else:
                    for step_rxn in new_rxn:
                        fout.write('{}\n'.format(step_rxn))
            except Exception as e:
                print(f'{i}th reaction has a problem of ...')
                print(e, '\n')
                pass

def generate_mechdata_unknown_condition(args):
    '''
    Input format is a string of 'reaction_smiles NameRXN_name'
    '''

    # Make directory for saving file
    save_dir_path=os.path.dirname(args.save)
    os.makedirs(save_dir_path, exist_ok=True)

    with open(args.data, 'r') as file, open(args.save, 'w') as fout:
        lines = file.readlines()
        for i, line in tqdm(enumerate(lines)):
            try:
                rxn = line.split()[0]
                if len(line.split()[1:])==1:
                    label = line.split()[1:]
                else:
                    label = ' '.join(line.split()[1:])
                rxn_dict = {'reaction_name': label,
                            'reaction_smiles': rxn,}
                rxn_dict = reagent_matching_for_single_reaction(rxn_dict, label)
                if rxn_dict['conditions']:
                    new_rxn = generate_mechanism_for_one_reaction(rxn_dict, args)
                    if args.all_info:
                        fout.write('{}\n'.format(new_rxn))
                    else:
                        for step_rxn in new_rxn:
                            fout.write('{}\n'.format(step_rxn))
            except Exception as e:
                if args.verbosity > 0:
                    print(f'{i}th reaction has a problem of ...')
                    print(e, '\n')
                pass