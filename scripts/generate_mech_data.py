import os
import logging
import json
import pickle
from tqdm import tqdm
from multiprocessing import Pool
from templates import Reaction_templates
from utils import reaction_process
from utils.elementary_reactions import Get_Reactions
from utils.exceptions import *
from utils.reaction_process import flatten_list
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

def merge_dicts(d, u):
    for k, v in u.items():
        if isinstance(v, dict):
            d[k] = merge_dicts(d.get(k, {}), v)
        else:
            d[k] = d.get(k, 0) + v
    return d

def calling_rxn_template(reaction_dict):
    """
    Retrieve the elementary reaction templates of each reaction
    reaction_dict: one reaction in reactions_with_conditions
    """
    rxn_class_name = reaction_dict['reaction_name']
    rxn_condition = reaction_dict['conditions']
    class_key = get_class_key(rxn_class_name)
    rxn_templates = Reaction_templates.class_reaction_templates[class_key]

    conditions = [rxn_templates[condition] for condition in rxn_condition]

    return conditions

def get_class_key(class_name):
    """
    returns the key(tuple) of class_reaction_templates that includes the class_name
    """

    for classes in Reaction_templates.class_reaction_templates.keys():
        if class_name in classes:
            return classes

    return None

def reagent_matching_for_single_reaction(reaction, class_key):
    '''
    To match reagents for a certain class from list of reaction, in the form of {'rxnsmiles': {'reaction_name': name}}.
    reactions: list of reactions
    '''
    reactants, agents, products = reaction['reaction_smiles'].split(">")
    mols = [Chem.MolFromSmiles(smi) for smi in reactants.split(".") + agents.split(".")]
    reaction['conditions'] = []
    class_key = get_class_key(class_key)
    if class_key in Reaction_templates.class_reaction_templates:
        class_data = Reaction_templates.class_reaction_templates[class_key]

        for cond_name, cond_data in class_data.items():
            if cond_data['Reagent'] or cond_data['Exclude_reagent']:
                cond_mols = []
                if cond_data['Reagent']:
                    cond_mols = [Chem.MolFromSmarts(x) for x in cond_data['Reagent']]
                exclude_cond_mols = []
                if cond_data['Exclude_reagent']:
                    exclude_cond_mols = [Chem.MolFromSmarts(x) for x in cond_data['Exclude_reagent']]

                matched_reagents = []
                exclude = False

                for patt in cond_mols:
                    for mol in mols:
                        if mol.GetSubstructMatch(patt):
                            matched_reagents.append(mol)
                            break
                for patt in exclude_cond_mols:
                    for mol in mols:
                        if mol.GetSubstructMatch(patt):
                            exclude = True
                            break
                    if exclude: break
                if not exclude and len(cond_mols) == len(matched_reagents):
                    reaction['conditions'].append(cond_name)
            else:
                reaction['conditions'].append(cond_name)
    return reaction

def get_mechanistic_network(rxn, args):
    conditions = calling_rxn_template(rxn)

    reaction_dict = {}
    statistic_dict = {}
    
    if args.verbosity:
        logging.info('Generating mechanism for: ' + rxn['reaction_smiles'])

    for i, cond in enumerate(conditions):
        reaction = reaction_process.Reaction_Network(rxn, i, args)
        reaction.set_template_dict(cond)
        first_error = 0
        for step in range(reaction.length):
            try:
                reaction.run_reaction(step)
            except NoReactantError as e:
                if first_error == 0:
                    first_error = 1
                continue
            except NoAcidBaseError as e:
                if first_error == 0:
                    first_error = 2
                continue
            except Exception as e:
                if first_error == 0:
                    first_error = 3
                if args.verbosity:
                    logging.info(f'{e}')
                continue

        if reaction.is_product_formed():
            elem_reactions = Get_Reactions(reaction)

            if not args.do_not_pruning:
                try:
                    elem_reactions.graph_pruning()
                except Exception as e:
                    statistic_dict[reaction.reaction_condition] = {'Reaction network pruning error': 1}
                    if args.verbosity:
                        logging.info(f'{e}')
                    continue
                try:
                    elem_reactions.get_elementary_reactions_info()
                except Exception as e:
                    statistic_dict[reaction.reaction_condition] = {'Getting reaction info. error': 1}
                    if args.verbosity:
                        logging.info(f'{e}')
                    continue

            try:
                rxn_smi = elem_reactions.convert_to_smiles()
            except Exception as e:
                statistic_dict[reaction.reaction_condition] = {'Conversion to SMILES error': 1}
                if args.verbosity:
                    logging.info(f'{e}')
                continue

            if args.all_info:
                reaction_dict[reaction.reaction_condition] = elem_reactions
                # print(elem_reactions.reaction_info)
            else:
                reaction_dict[reaction.reaction_condition] = rxn_smi
                
            statistic_dict[reaction.reaction_condition] = {'Success': {'overall reactions': 1,
                                        'elementary reactions': len(rxn_smi),
                                        'reactants': len(elem_reactions.reactant_node),
                                        'products': len(elem_reactions.product_node),
                                        'byproducts': len(elem_reactions.byproduct_node),
                                        'intermediates': len(elem_reactions.intermediate_node),
                                        'spectators': len(elem_reactions.spectator_node),
                                        }}
            if args.verbosity: 
                logging.info(f'Product is formed for {reaction.reaction_condition}, total {len(rxn_smi)} elem. steps were generated ')

                
            # elem_reactions.print_graph()


        else:
            if first_error == 0:
                statistic_dict[reaction.reaction_condition] = {'No product': 1}
                if args.verbosity: 
                    logging.info(f'Product is not produced for {reaction.reaction_condition}')
            elif first_error == 1:
                statistic_dict[reaction.reaction_condition] = {'No reactant': 1}
                if args.verbosity: 
                    logging.info(f'Product is not produced for {reaction.reaction_condition} due to no reactant')
            elif first_error == 2:
                statistic_dict[reaction.reaction_condition] = {'No acid base': 1}
                if args.verbosity: 
                    logging.info(f'Product is not produced for {reaction.reaction_condition} due to no acids or bases')
            elif first_error == 3:
                statistic_dict[reaction.reaction_condition] = {'Error': 1}
                if args.verbosity: 
                    logging.info(f'Error occured for {reaction.reaction_condition}')

    return reaction_dict, statistic_dict

def generate_elementary_reaction(input):
    """
    It makes elementary reactions for a single overall reactions.
    It takes a line of 'reactants>>products NameRXN_class'
    It returns two dictionaries for the reaction networks and statistics
    """
    # try:
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

    try:
        rxn_dict = reagent_matching_for_single_reaction(rxn_dict, label)
        if not rxn_dict['conditions']:
            statistics = {'No templates' : 1}
            return rxn, [], {label: statistics}

        results, statistics = get_mechanistic_network(rxn_dict, args)
    except Exception as e:
        if args.verbosity:
            logging.info(f'{e}')
        statistics = {'Invalid reaction' : 1}
        return rxn, [], {label: statistics}

    return rxn, list(results.values()), {label: statistics}
    
def find_elementary_reactions(d):
    if isinstance(d, dict):
        for key, value in d.items():
            if key == 'elementary reactions':
                return value
            elif isinstance(value, dict):
                result = find_elementary_reactions(value)
                if result is not None:
                    return result
    return None

def generate_mechdata_multiprocess(args):
    '''
    Input format is a string of 'reaction_smiles NameRXN_name'
    '''

    # Make directory for saving file
    save_dir_path=os.path.dirname(args.save)
    os.makedirs(save_dir_path, exist_ok=True)

    p = Pool(args.process)
    num_generated_rxn=0
    num_used_rxn = 0

    with open(args.data, 'r') as file:
        lines = file.readlines()
        iterables = [(line, args) for line in lines]

    if args.stat:
        statistics = {}


    base_file_root, _ = os.path.splitext(args.save)
    debug_file_path = f"{base_file_root}_debug.txt"

    if args.all_info:
        writing_format = 'wb'
    else: writing_format = 'w'

    with open(args.save, writing_format) as fout, open(debug_file_path, 'w') as file_debug:
        for results in tqdm(p.imap_unordered(generate_elementary_reaction, iterables), total=len(lines)):
            if results is not None:
                rxn, elem_rxns, stat = results
                if elem_rxns:
                    num_used_rxn += 1
                    if args.all_info:
                        if find_elementary_reactions(stat):
                            num_generated_rxn += find_elementary_reactions(stat)
                        pickle.dump(elem_rxns, fout, protocol=pickle.HIGHEST_PROTOCOL)
                    else:
                        elem_rxns = flatten_list(elem_rxns)
                        num_generated_rxn += len(elem_rxns)
                        fout.write('\n'.join(elem_rxns) + '\n')

                else:
                    file_debug.write(f"{rxn} {stat}\n")
                    
                if args.stat:
                    merge_dicts(statistics, stat)

        if args.stat:
            base_file_root, _ = os.path.splitext(args.save)
            stat_file_path = f"{base_file_root}_statistics.json"
            with open(stat_file_path, 'w') as file:
                json.dump(statistics, file, indent=4)


    if not args.debug:
        os.remove(debug_file_path)

    logging.info(f'In total, {num_used_rxn} overall reactions were utilized')
    logging.info(f'{num_generated_rxn} elementary steps were generated')
    logging.info(f'Mechanistic dataset generation is done')

