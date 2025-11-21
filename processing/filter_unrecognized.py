import random

# in_file = "../data/unrecognized.txt"
# out_file = "../data/unrecognized_filtered.txt"

# with open (in_file, 'r') as f, open (out_file, 'w') as o:
#     for line in f:
#         if "Pd" in line or "Cu" in line:
#             o.write(line)

# in_file = "../data/unrec_f_p_d.txt"
# out_file = "../data/unrec_one_p.txt"

# with open (in_file, 'r') as f, open (out_file, 'w') as o:
#     for line in f:
#         parts = line.split(">>")
#         if len(parts[1].split(".")) == 1:
#             o.write(line.split("|")[0] + "\n")

# in_file = "../data/unrec_one_p.txt"
# out_file = "../data/unrec_5_percent.txt"

# seed = 42  # 👈 change this to any integer you like

# # ======================
# random.seed(seed)

# # Read all lines
# with open(in_file, 'r', encoding='utf-8') as f:
#     lines = f.readlines()

# # Compute 1% of total lines (at least one)
# subset_size = max(1, len(lines) // 100) * 5

# # Randomly sample 1% of lines
# subset = random.sample(lines, subset_size)

# # Write subset to new file
# with open(out_file, 'w', encoding='utf-8') as f:
#     f.writelines(subset)

in_file = "../data/unrec_5_percent.txt"
out_file = "../data/unrec_5_trimmed.txt"

# ======================

# Read all lines
with open(in_file, 'r', encoding='utf-8') as f:
    lines = f.read().splitlines()

# Sort lines by length (longest first)
sorted_lines = sorted(lines, key=len, reverse=True)

# Determine how many to drop (top 50%)
n = len(sorted_lines)
print(n)

lower_cut = n // 4
upper_cut = n - lower_cut
print(lower_cut, upper_cut)

kept_lines = sorted_lines[lower_cut:upper_cut]

# Write them back, preserving one per line
with open(out_file, 'w', encoding='utf-8') as f:
    for line in kept_lines:
        f.write(line + '\n')