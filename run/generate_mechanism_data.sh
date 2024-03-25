python generate_mechanism_data.py \
  --data 'data/uspto_classified.txt'\
  --save 'results/uspto_all.jsonl'\
  --proton True\
  --byproduct False\
  --spectator False\
  --end False\
  --debug True\
  --verbosity 0\
  --process 30\
  --all_info True\
  #--do_not_pruning True
  #  --data 'data/uspto_classified.txt'\
#  --rxn_class '"CO2H-Et deprotection":'\