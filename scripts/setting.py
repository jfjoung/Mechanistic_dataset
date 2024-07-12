class Args:
    def __init__(self):
        self.num_combination = 30
        self.uni_rxn = True
        self.proton = True
        self.max_num_temp = 50
        self.stoichiometry = True
        self.simple = False
        self.do_not_pruning = False
        self.num_cycles = 9
        self.num_reaction_node = 50
        self.byproduct = True
        self.spectator = True
        self.full = False
        self.end = True
        self.plain = False
        self.explicit_H = True
        self.reagent = False
        self.remapping = True
        self.data = './data/test_data.txt'
        self.save = './results2/test.txt'
        self.debug = False
        self.all_info = False
        self.stat = True
        self.rxn_class = None
        self.process = 10
        self.verbosity = 0
  