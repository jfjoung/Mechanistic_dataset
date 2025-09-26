from templates.Reaction_templates import class_reaction_templates

all_templates = []
all_a_b = []

for val in class_reaction_templates.values():
    for rxn in val:
        stages = rxn.value()['Stages']
        ###get all templates


###process all templates w/ a/b

### a/b process, only look in current reactant dict

for batch in reactions: 
    randomize templates order
    try template against reactants 
    then, write out products, maybe do 2 templates per reaction?


worker process:
    takes batch, processes and writes


for lines in reactions:
    process a batch and send to worker