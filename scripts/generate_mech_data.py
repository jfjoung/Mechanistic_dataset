import os
import sys
import logging
import networkx as nx
from tqdm import tqdm
import json
from utils.apply_mechanistic_template import get_mechanistic_network, elementary_reaction, flatten_list, reagent_matching_for_single_reaction


def generate_mechanism_for_one_reaction(rxn, args):
    if args.verbosity > 0:
        logging.info('Generating mechanism for: ' + rxn['reaction_smiles'])
    G_dict = get_mechanistic_network(rxn, args)
    elem_steps_stats = {} # To get statistics
    if args.all_info:
        elem_dict = dict()
    else:
        elem_list = list()

    for cond, G in G_dict.items():
        try:
            elem_rxns = elementary_reaction(G, args)
            elem_steps_stats[cond] = len(elem_rxns)
            if args.all_info:
                elem_dict[cond]={'Reaction graph': nx.node_link_data(G),
                                 'Elementary steps': elem_rxns}
            else:
                elem_list.append(elem_rxns)

        except Exception as e:
            if args.verbosity > 0:
                logging.info('Error occured in {}'.format(cond))
                logging.info(f'{e}')
            pass

    if args.all_info:
        rxn['Mechanism'] = elem_dict
        return rxn, elem_steps_stats
    else:
        return flatten_list(elem_list), elem_steps_stats

def generate_mechdata(args):
    '''
    Input format is a string of 'reaction_smiles NameRXN_name'
    '''

    # Make directory for saving file
    save_dir_path=os.path.dirname(args.save)
    os.makedirs(save_dir_path, exist_ok=True)

    conditions_stats = {}  # Statistics on how many reaction is met the reaction condition
    elem_steps_per_cond = {} # Statistics on how many elementary reaction is generated
    num_generated_rxn=0

    with open(args.data, 'r') as file, open(args.save, 'w') as fout:
        lines = file.readlines()
        for i, line in enumerate(tqdm(lines)):
            if args.verbosity > 0:
                logging.info(f'Started {i}th reaction')

            rxn = line.split()[0]
            if len(line.split()[1:])==1:
                label = line.split()[-1]
            else:
                label = ' '.join(line.split()[1:])
            rxn_dict = {'reaction_name': label,
                        'reaction_smiles': rxn,}
            if label not in conditions_stats:
                conditions_stats[label] = {'RXN_generated': 0, 'Failed': 0,
                                           'no_template': 0}

            rxn_dict = reagent_matching_for_single_reaction(rxn_dict, label)
            if rxn_dict['conditions']:
                try:
                    new_rxn, elem_steps_stats = generate_mechanism_for_one_reaction(rxn_dict, args)
                    if new_rxn:
                        conditions_stats[label]['RXN_generated'] += 1
                        if label not in elem_steps_per_cond:
                            elem_steps_per_cond[label] = {}
                        for cond, steps in elem_steps_stats.items():
                            if cond in elem_steps_per_cond[label]:
                                elem_steps_per_cond[label][cond] += steps
                            else:
                                elem_steps_per_cond[label][cond] = steps
                    else:
                        conditions_stats[label]['Failed'] += 1
                    if args.all_info:
                        fout.write('{}\n'.format(new_rxn))
                    else:
                        if args.verbosity > 0:
                            logging.info(
                                'Total {} elementary reactions are generated \n'.format(len(new_rxn)))
                        for step_rxn in new_rxn:
                            num_generated_rxn += 1
                            fout.write('{}\n'.format(step_rxn))

                except Exception as e:
                    conditions_stats[label]['Failed'] += 1
                    if args.verbosity > 0:
                        logging.info('The reaction has problem!')
                        logging.info(f'{line}')
                        logging.info(f'{e}\n')
                    pass
            else:
                conditions_stats[label]['no_template'] += 1

    if args.stat:
        base_file_root, _ = os.path.splitext(args.save)
        condition_stat_file_path = f"{base_file_root}_condition_stat.json"
        elementary_stat_file_path = f"{base_file_root}_elementary_stat.json"
        with open(condition_stat_file_path, 'w') as file:
            json.dump(conditions_stats, file, indent=4)
        with open(elementary_stat_file_path, 'w') as file:
            json.dump(elem_steps_per_cond, file, indent=4)

    logging.info(f'{num_generated_rxn} reactions generated')
    logging.info('The generation process has ended')