python generate_mechanism_data.py \
  --data 'data/USPTO_train.txt'\
  --save 'results/USPTO_train.txt'\
  --verbosity 0\
  --proton True\
  --byproduct True\
  --spectator True\
  --reagent False\
  --full False\
  --end True\
  --debug False\
  --explicit_H True\
  --process 30\
  --all_info False\


python generate_mechanism_data.py \
  --data 'data/USPTO_valid.txt'\
  --save 'results/USPTO_valid.txt'\
  --verbosity 0\
  --proton True\
  --byproduct True\
  --spectator True\
  --reagent False\
  --full False\
  --end True\
  --debug False\
  --explicit_H True\
  --process 30\
  --all_info False\

python generate_mechanism_data.py \
  --data 'data/USPTO_test.txt'\
  --save 'results/USPTO_test.txt'\
  --verbosity 0\
  --proton True\
  --byproduct True\
  --spectator True\
  --reagent False\
  --full False\
  --end True\
  --debug False\
  --explicit_H True\
  --process 30\
  --all_info False\

  #  --remapping True\
  #--do_not_pruning True\

#  --data 'data/testing.txt'\
#  --save 'results/testing.jsonl'\
#  --verbosity 4\
  #  --data 'data/uspto_classified.txt'\
#  --rxn_class '"CO2H-Et deprotection":'\
