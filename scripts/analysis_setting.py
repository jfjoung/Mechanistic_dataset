class Args:
    def __init__(self):
        self.data = './results/USPTO.pickle'  # Path to the file where results were saved.
        self.template_analysis = True
        self.ERS_analysis = True
        self.mol_analysis = True
        self.process = 32  # Number of processors to use in multiprocessing.
