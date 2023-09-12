# ======================================================================================================================
# Yingying Dong. 2023-9. Getting the codon usage frequency of every specie.
# Details: Removing useless columns and splitting the large table (all species CUB) into multiple files.
# =====================================================================================================================
import os

if os.path.exists('./all_species'):
    print('exists directory ')
else:
    os.mkdir('./all_species')
input_file = 'o537-Refseq_species.tsv'
output_file = 'refseq_species.txt'

columns_to_delete = list(range(1, 4)) + list(range(5, 13))

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        columns = line.strip().split('\t')

        # Create a new line with the specified columns removed
        new_line = [columns[i - 1] for i in range(1, len(columns) + 1) if i not in columns_to_delete]

        # Write the new line to the output file
        outfile.write("\t".join(new_line) + "\n")

with open(output_file, 'r') as all_spe:
    # Read the header line
    header = all_spe.readline().strip().split('\t')
    # Initialize a dictionary to store file handles
    file_handles = {}
    # Read the rest of the lines
    for line in all_spe:
        # Split the line into fields
        fields = line.strip().split('\t')

        # Get the first column (ID)
        id_value = fields[0]
        print(id_value)

        '''# Check if a file handle already exists for this ID
        if id_value in file_handles:
            # If it exists, write the line to the existing file
            file_handles[id_value].write(line)
        else:
        # If it doesn't exist, create a new file handle and write the header and line'''

        file_handles[id_value] = open(f'./all_species/{id_value}.txt', 'w')
        file_handles[id_value].write('\t'.join(header) + '\n')
        file_handles[id_value].write(line)
        file_handles[id_value].close()

'''# Close all file handles
for file_handle in file_handles.values():
    file_handle.close()'''
