import os
import sys
import logging
from tqdm import tqdm
from multiprocessing import Pool
from utils import Reaction_templates, reaction_process
from utils.elementary_reactions import Get_Reactions
from utils.exceptions import *
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')



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
    first_error = 0

    for i, cond in enumerate(conditions):
        reaction = reaction_process.Reaction_Network(rxn, i, args)
        reaction.set_template_dict(cond)
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
                print(e)
                if first_error == 0:
                    first_error = 3
                continue

        if reaction.is_product_formed():
            elem_reactions = Get_Reactions(reaction)

            if not args.do_not_pruning:
                try:
                    elem_reactions.graph_pruning()
                except:
                    statistic_dict[reaction.reaction_condition] = {'Reaction network pruning error': 1}
                    print(rxn)
                    raise
                    continue
            elem_reactions.get_elementary_reactions_info()

            if args.all_info:
                reaction_dict[reaction.reaction_condition] = elem_reactions
            else:
                rxn_smi = elem_reactions.convert_to_smiles()
                reaction_dict[reaction.reaction_condition] = rxn_smi
            # elem_reactions.print_graph()
            statistic_dict[reaction.reaction_condition] = {'overall reactions': 1,
                                                          'elementary reactions': len(elem_reactions.reaction_node),
                                                           'reactants': len(elem_reactions.reactant_node),
                                                           'products': len(elem_reactions.product_node),
                                                           'byproducts': len(elem_reactions.byproduct_node),
                                                           'intermediates': len(elem_reactions.intermediate_node),
                                                           'spectators': len(elem_reactions.spectator_node),
                                                           }
        else:
            if first_error == 0:
                statistic_dict[reaction.reaction_condition] = {'No product': 1}
            elif first_error == 1:
                statistic_dict[reaction.reaction_condition] = {'No reactant': 1}
            elif first_error == 2:
                statistic_dict[reaction.reaction_condition] = {'No acid base': 1}
            elif first_error == 3:
                statistic_dict[reaction.reaction_condition] = {'Error': 1}

    return reaction_dict, statistic_dict

def generate_elementary_reaction(input):
    """
    It makes elementary reactions for a single overall reactions.
    It takes a line of 'reactants>>products NameRXN_class'
    It returns two dictionaries for the reaction networks and statistics
    """
    line, args = input
    reaction_networks = {}


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

    rxn_dict = reagent_matching_for_single_reaction(rxn_dict, label)
    if not rxn_dict['conditions']:
        statistics = {'No templates' : 1}
        return rxn, {label: reaction_networks}, {label: statistics}

    results, statistics = get_mechanistic_network(rxn_dict, args)

    return rxn, {label: results}, {label: statistics}



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

    if args.stat:
        statistics = {}


    base_file_root, _ = os.path.splitext(args.save)
    debug_file_path = f"{base_file_root}_debug.txt"


    with open(args.save, 'w') as fout, open(debug_file_path, 'w') as file_debug:
        iterables = [(line, args) for line in lines]
        for result in tqdm(p.imap_unordered(generate_elementary_reaction, iterables), total=len(lines)):
            rxn, results, stat = result
            print(rxn)
            print(stat, '\n')

            # if args.debug:
            #     failed_reaction = []
            # if not result: continue
            # if result[0]:
            #     if args.debug:
            #         rxn, results, stat = result
            #     else:
            #         _, elem_reaction, stat = result
            #     if elem_reaction:
            #         num_used_rxn += 1
            #         if args.all_info:
            #             fout.write(json.dumps(elem_reaction) + "\n")
            #             # json.dump(elem_reaction, fout, indent=4)
            #             # fout.write("\n")
            #         else:
            #             fout.write('\n'.join(elem_reaction) + '\n')
            #     num_generated_rxn += len(elem_reaction)
            #     if args.stat:
            #         merge_dicts(statistics, stat)
            #     if args.debug:
            #         for failed in failed_products:
            #             failed_reaction.append(';'.join(failed)+ '\n')
            # else:
            #     _, error, line, stat = result
            #     if args.stat:
            #         merge_dicts(statistics, stat)
            #     if args.debug:
            #         failed_reaction.append(f'{line};No condition;{error}'+ '\n')
            # if args.debug:
            #     if failed_reaction:
            #         file_debug.write(failed_reaction[0])
