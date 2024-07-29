import pandas as pd
import numpy as np

# Load the file with newline as delimiter
file_path = 'uspto_classified.txt'
with open(file_path, 'r') as file:
    data = file.read().split('\n')

# Convert the data to a DataFrame
data = pd.DataFrame(data, columns=['text'])

# Define a random seed for reproducibility
random_seed = 1205
np.random.seed(random_seed)

# Shuffle the data
shuffled_data = data.sample(frac=1, random_state=random_seed).reset_index(drop=True)

# Calculate the split indices
train_end = int(0.8 * len(shuffled_data))
valid_end = int(0.9 * len(shuffled_data))

# Split the data
train_data = shuffled_data[:train_end]
valid_data = shuffled_data[train_end:valid_end]
test_data = shuffled_data[valid_end:]

# Function to save data to a text file
def save_to_textfile(dataframe, filename):
    with open(filename, 'w') as f:
        for text in dataframe['text']:
            f.write(f"{text}\n")

# Save the datasets into separate text files
save_to_textfile(train_data, 'USPTO_train.txt')
save_to_textfile(valid_data, 'USPTO_valid.txt')
save_to_textfile(test_data, 'USPTO_test.txt')

# Display the first few rows of each dataset
print("Train Data:")
print(train_data.head())
print("\nValidation Data:")
print(valid_data.head())
print("\nTest Data:")
print(test_data.head())