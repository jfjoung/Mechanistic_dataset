import random

input_file = "../finalized_data/test_d.txt"
output_file = "../finalized_data/test_d_sub.txt"
sample_fraction = 0.05  # 5%
random_seed = 42  # change this for different reproducible results

random.seed(random_seed)

batches = {}
standalone_strings = []

# Step 1: Read and group
with open(input_file, "r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        parts = line.split("|", 1)
        if len(parts) == 2 and parts[1].isdigit():  # |number
            num = parts[1]
            batches.setdefault(num, []).append(line)
        else:
            standalone_strings.append(line)

# Step 2: Create unified list of "units"
units = list(batches.values()) + [[s] for s in standalone_strings]

# Step 3: Randomly pick 5%
sample_size = max(1, int(len(units) * sample_fraction))
sampled_units = random.sample(units, sample_size)

# Step 4: Flatten and write output
with open(output_file, "w", encoding="utf-8") as f:
    for unit in sampled_units:
        for line in unit:
            f.write(line + "\n")
