import os
import sys
import logging
import networkx as nx
from tqdm import tqdm
import json
from multiprocessing import Pool
from utils.apply_mechanistic_template import get_mechanistic_network, elementary_reaction, flatten_list, reagent_matching_for_single_reaction, find_chemical_nodes

def generate_mechdata_single(input):
    line, args = input

    rxn = line.split()[0]
    if len(line.split()[1:]) == 1:
        label = line.split()[-1]
    else:
        label = ' '.join(line.split()[1:])

    # When reaction class is specified, then corresponding reactions will only be considered.
    if args.rxn_class and args.rxn_class not in label:
        return None

    rxn_dict = {'reaction_name': label,
                'reaction_smiles': rxn, }

    elem_steps_stats = {label: {'No templates': 0,
                                'Error': 0,}}  # To get statistics

    rxn_dict = reagent_matching_for_single_reaction(rxn_dict, label)
    if not rxn_dict['conditions']:
        elem_steps_stats[label]['No templates']+=1
        return elem_steps_stats
    else:
        if args.verbosity > 0:
            logging.info('Generating mechanism for: ' + rxn_dict['reaction_smiles'])
        try:
            G_dict = get_mechanistic_network(rxn_dict, args)
        except Exception as e:
            elem_steps_stats[label]['Error'] += 1
            if args.verbosity > 0:
                logging.info('Error occured in {}'.format(rxn))
                logging.info(f'{e}')
            return elem_steps_stats  # 0

        if args.all_info:
            elem_dict = dict()
        else:
            elem_list = list()

        for cond, G in G_dict.items():
            elem_steps_stats[label][cond]= {'No products': 0,
                                            'Too many reactions': 0,
                                            'Error': 0,
                                            'Success': {}}

            if G == "Products are not produced.":
                elem_steps_stats[label][cond]['No products'] += 1
                if args.verbosity > 0:
                    logging.info('Products are not produced for the reaction of {}'.format(cond))

            else:
                try:
                    elem_rxns = elementary_reaction(G, args)
                    reaction_node, reactant_node, product_node, byproduct_node, intermediate_node, spectator_node=find_chemical_nodes(G)
                    elem_steps_stats[label][cond]['Success'] = {'Num of overall reaction': 1,
                                                                 'Num of total elem steps': len(elem_rxns),
                                                                 'Num of total reactants': len(reactant_node),
                                                                 'Num of total products': len(product_node),
                                                                 'Num of total byproducts': len(byproduct_node),
                                                                 'Num of total intermediates': len(intermediate_node),
                                                                 'Num of total spectators': len(spectator_node),
                                                                 }

                    if args.all_info:
                        elem_dict[cond] = {'Reaction graph': nx.node_link_data(G),
                                           'Elementary steps': elem_rxns}
                    else:
                        elem_list.append(elem_rxns)

                except Exception as e:
                    if str(e) == "Too many reaction nodes.":
                        if args.verbosity > 0:
                            logging.info('There are too many steps for the reaction of {}'.format(cond))
                        elem_steps_stats[label][cond]['Too many reactions']+=1
                        # return elem_steps_stats
                    elif str(e) == "Too many cycles.":
                        if args.verbosity > 0:
                            logging.info('There are too catalytic cycles for the reaction of {}'.format(cond))
                        elem_steps_stats[label][cond]['Too many reactions']+=1
                        # return elem_steps_stats
                    else:
                        if args.verbosity > 0:
                            logging.info('Error occured in {}'.format(cond))
                            logging.info(f'{e}')
                        elem_steps_stats[label][cond]['Error'] += 1
                        # return elem_steps_stats

        if args.all_info and elem_dict:
            rxn_dict['Mechanism'] = elem_dict
            return rxn_dict, elem_steps_stats
        elif elem_list:
            return flatten_list(elem_list), elem_steps_stats
        else:
            return elem_steps_stats

def merge_dicts(d, u):
    for k, v in u.items():
        if isinstance(v, dict):
            d[k] = merge_dicts(d.get(k, {}), v)
        else:
            d[k] = d.get(k, 0) + v
    return d

def generate_mechdata_multiprocess(args):
    '''
    Input format is a string of 'reaction_smiles NameRXN_name'
    '''

    # Make directory for saving file
    save_dir_path=os.path.dirname(args.save)
    os.makedirs(save_dir_path, exist_ok=True)

    p = Pool(args.process)
    num_generated_rxn=0

    with open(args.data, 'r') as file:
        lines = file.readlines()

    if args.stat:
        statistics = {}

    with open(args.save, 'w') as fout:
        iterables = [(line, args) for line in lines]
        for result in tqdm(p.imap(generate_mechdata_single, iterables), total=len(lines)):
            if isinstance(result, dict):
                if args.stat:
                    merge_dicts(statistics, result)
            elif result:
                elem_reaction, stat = result
                fout.write('\n'.join(elem_reaction) + '\n')
                num_generated_rxn+=len(elem_reaction)
                if args.stat:
                    merge_dicts(statistics, stat)

    if args.stat:
        base_file_root, _ = os.path.splitext(args.save)
        stat_file_path = f"{base_file_root}_statistics.json"
        with open(stat_file_path, 'w') as file:
            json.dump(statistics, file, indent=4)


    logging.info(f'{num_generated_rxn} reactions generated')
    logging.info('The generation process has ended')
