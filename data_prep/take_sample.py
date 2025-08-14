import random

input_file = '../test_runs/full_suzuki_prepped.txt'
output_file = '../test_runs/suzuki_sample.txt'
sample_size = 10

# Step 1: Read all lines
with open(input_file, 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Step 2: Sample 10,000 lines without replacement
if len(lines) < sample_size:
    raise ValueError(f"File only has {len(lines)} lines, cannot sample {sample_size} unique lines.")

sampled_lines = random.sample(lines, sample_size)

# Step 3: Write sampled lines to output file
with open(output_file, 'w', encoding='utf-8') as f:
    f.writelines(sampled_lines)

print(f"Sampled {sample_size} lines written to {output_file}")