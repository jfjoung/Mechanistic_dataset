class Args:
    def __init__(self):
        self.num_combination = 30  # Limits the number of possible reactant combinations per template to prevent combinatorial explosion.
        self.uni_rxn = True  # Adds a unimolecular version of a bimolecular reaction template while retaining the original bimolecular template.
        self.proton = True  # Ensures proton balance; if True, the program checks if an acid or base from AcidBase_Lookup is present in the reaction flask and adjusts the template on-the-fly to maintain proton balance.
        self.max_num_temp = 30  # Limits the number of possible templates when uni_rxn or proton is set to True, preventing combinatorial explosion.
        self.stoichiometry = True  # Allows adding duplicated reactants to the flask if the stoichiometry of a reactant used more than once is recorded incorrectly.
        self.do_not_pruning = False  # Determines whether to prune the reaction network, leaving only those that produce the final product.
        self.num_cycles = 9  # Sets an upper limit on the number of catalytic cycles in the reaction network to prevent excessive computation time when converting to elementary reaction SMARTS.
        self.num_reaction_node = 50  # Sets an upper limit on the number of reactions in the reaction network to prevent excessive computation time when converting to elementary reaction SMARTS.
        self.byproduct = True  # Adds byproducts generated in previous reactions when converting to elementary reaction SMARTS.
        self.spectator = True  # Adds all chemical species present in the flask at the current step when converting to elementary reaction SMARTS.
        self.full = False  # Determines whether to return overall reaction SMARTS instead of elementary reaction SMARTS.
        self.end = True  # Adds reactions that generate a product from an already existing product.
        self.plain = False  # Generates reaction SMARTS without atom mapping numbers.
        self.explicit_H = False  # Explicitly represents hydrogens in the reaction SMARTS.
        self.reagent = True  # Places non-participating chemical species between '>>' in Reaction SMARTS formatted as Reactants > Reagents > Products.
        self.remapping = True  # Remaps the atom mapping in Reaction SMARTS to start from 1.
        self.data = './data/uspto_classified.txt'  # Path to the file containing overall reactions.
        self.save = './results/USPTO.pickle'  # Path to the file where results will be saved.
        self.debug = True  # Records reactions with errors for debugging purposes.
        self.all_info = True  # Saves all reactions, including the reaction network.
        self.stat = True  # Generates statistics for the reactions.
        self.rxn_class = None  # Specifies a particular reaction class to extract elementary reactions from.
        self.process = 32  # Number of processors to use in multiprocessing.
        self.verbosity = False  # Sets the logging.