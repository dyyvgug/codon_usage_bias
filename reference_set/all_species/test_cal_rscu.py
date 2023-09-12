# Open the file for reading
with open('Aaosphaeria arxii CBS 175.79.txt', 'r') as file:
    # Read the lines from the file
    lines = file.readlines()

# Split the lines into two rows

table = [line.strip().split('\t') for line in lines]
for row in table:
    del row[0]

# Extract keys and values
keys = table[0]
values = table[1]

# Create a dictionary using a dictionary comprehension
result_dict = {key: value for key, value in zip(keys, values)}

print(result_dict)
