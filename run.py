import os
import logging
from scripts.setting import Args
from datetime import datetime
from scripts.generate_mech_data import generate_mechdata_multiprocess


if __name__ == '__main__':
    '''
    You need a reaction string of 'reaction_smiles NameRXN_name', 
    generate_mechdata_multiprocess(args)
    '''

    args=Args()
    os.makedirs('./logs/', exist_ok=True)
    time = datetime.now().strftime('%Y-%m-%d_%H-%M')
    logging.basicConfig(filename=f'./logs/generate_mechanism_data_{time}.log', level=logging.DEBUG, format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.info('Starting mechanistic dataset generation...')
    logging.info(f'Arguments')
    for key, value in vars(args).items():
        logging.info('{}: {}'.format(key, value))

    generate_mechdata_multiprocess(args)

    # args.data = './data/uspto_val.txt'
    # args.save = './results/flow_ood/Negishi_coupling_val.txt'

    # generate_mechdata_multiprocess(args)

    # args.data = './data/uspto_test.txt'
    # args.save = './results/flow_ood/Negishi_coupling_test.txt'

    # generate_mechdata_multiprocess(args)