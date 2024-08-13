class Args:
    def __init__(self):
        self.data = './results/test.pickle'  # Path to the file where results were saved.
        self.template_analysis = False
        self.ERS_analysis = False
        self.mol_analysis = False
        self.process = 32  # Number of processors to use in multiprocessing.
        self.validity = True
