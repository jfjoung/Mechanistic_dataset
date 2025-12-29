import numpy as np

# File paths
input_path = '../results/templated_d.txt'
train_path = '../data/barluenga/train.txt'
valid_path = '../data/barluenga/val.txt'
test_path  = '../data/barluenga/test.txt'

# Set seed and initialize
random_seed = 1205
np.random.seed(random_seed)

# Step 1: Get all unique numbers from the file
print("Indexing unique numbers...")
unique_numbers = set()
with open(input_path, 'r') as f:
    for line in f:
        line = line.strip()
        if '|' not in line:
            continue
        _, num = line.rsplit('|', 1)
        unique_numbers.add(num)

# Step 2: Shuffle and split numbers
unique_numbers = list(unique_numbers)
np.random.shuffle(unique_numbers)
n = len(unique_numbers)

train_keys = set(unique_numbers[:32])
valid_keys = set(unique_numbers[32:40])
test_keys  = set(unique_numbers[40:])

# Step 3: Stream through input file and write to appropriate output
print("Streaming and writing output...")
with open(input_path, 'r') as fin, \
     open(train_path, 'w') as ftrain, \
     open(valid_path, 'w') as fvalid, \
     open(test_path,  'w') as ftest:
    
    for line in fin:
        line = line.strip()
        if '|' not in line:
            continue
        _, num = line.rsplit('|', 1)

        if num in train_keys:
            ftrain.write(line + '\n')
        elif num in valid_keys:
            fvalid.write(line + '\n')
        elif num in test_keys:
            ftest.write(line + '\n')
