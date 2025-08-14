with open("../results/untemplated.txt", "r") as infile, open("../results/untemplated_fixed.txt", "w") as outfile:
    for line in infile:
        if line.strip():  # Only write non-blank lines
            outfile.write(line)