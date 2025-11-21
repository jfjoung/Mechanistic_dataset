import pandas as pd
import csv
import pyarrow.parquet as pq
import re

# Parameters
input_file = "../outputs/p_labeled.parquet"             # Path to input parquet file
output_file = "../data/o_allyl.txt"      # Output .txt file
allowed_values = {'O-allyl deprotection'} 

parquet_file = pq.ParquetFile(input_file)

def clean_text(content):
    """Remove all occurrences of |...| and any trailing space."""
    return re.sub(r'\|.*?\|\s*', '', content)

# Get column names
columns = parquet_file.schema.names
if len(columns) < 2:
    raise ValueError("The input file must have at least two columns.")
second_col = columns[1]

with open(output_file, "w") as f:
    for i in range(parquet_file.num_row_groups):
        # Read one row group at a time
        table = parquet_file.read_row_group(i)
        df = table.to_pandas()
        filtered_df = df[df[second_col].isin(allowed_values)]
        for row in filtered_df.itertuples(index=False):
            f.write(" ".join(clean_text(str(x)) for x in row) + "\n")
