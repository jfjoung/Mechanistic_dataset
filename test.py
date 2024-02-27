import os
import sys
import networkx as nx
from utils.parsing import parse_arguments
from scripts.generate_mech_data import generate_batch
import json


if __name__ == '__main__':
    args = parse_arguments()
    print(args)
    generate_batch(args)

