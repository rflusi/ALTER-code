import os

# Parameters
input_filename = "NoLinker_Mutant_MariaInput_V1.txt"  # input file
output_prefix = "NoLinker_V1"  # base name for output
max_samples_per_chunk = 500000  # maximum samples per file (adjust if you want)

# Load input
with open(input_filename, 'r') as f:
    lines = f.readlines()

# Remove header if there is one
if lines[0].startswith("peptide") or lines[0].startswith("Peptide"):
    header = lines[0]
    lines = lines[1:]
else:
    header = None

# Make split_inputs folder if needed
if not os.path.exists('split_NoLinker_V1'):
    os.makedirs('split_NoLinker_V1')

# Split into chunks
num_chunks = (len(lines) + max_samples_per_chunk - 1) // max_samples_per_chunk

for i in range(num_chunks):
    chunk = lines[i * max_samples_per_chunk : (i+1) * max_samples_per_chunk]
    chunk_filename = os.path.join('split_NoLinker_V1', '{}_chunk_{}.txt'.format(output_prefix, i+1))
    with open(chunk_filename, 'w') as out_f:
        if header:
            out_f.write(header)
        out_f.writelines(chunk)

print("Splitting completed. {} chunks created.".format(num_chunks))
