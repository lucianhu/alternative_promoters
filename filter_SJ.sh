#!/bin/bash

# Check if input and output folders are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 input_folder output_folder"
    exit 1
fi

input_folder=$1
output_folder=$2

# Check if input folder exists
if [ ! -d "$input_folder" ]; then
    echo "Input folder $input_folder not found!"
    exit 1
fi

# Create output folder if it doesn't exist
mkdir -p "$output_folder"

# Loop over all .SJ.out.tab files in the input folder
for file in "$input_folder"/*.SJ.out.tab; do
    # Check if there are any files in the folder
    if [ ! -f "$file" ]; then
        echo "No .SJ.out.tab files found in $input_folder!"
        exit 1
    fi

    # Extract the base name of the file (e.g., HIGH_96672 from HIGH_96672.SJ.out.tab)
    base_name=$(basename "$file" .SJ.out.tab)

    # Apply the Python filter script to the file and save the output to the output folder
    python3 - <<END
import pandas as pd

# Define columns
columns = ['chromosome', 'intron_start', 'intron_end', 'strand', 'intron_motif', 
           'annotated', 'unique_reads', 'multimapping_reads', 'max_overhang']

# Read the input file
file_name = '$file'
df = pd.read_csv(file_name, sep='\t', header=None, names=columns)

# Calculate intron length
df['intron_length'] = df['intron_end'] - df['intron_start']

# Define the filtering function
def filter_junctions(row):
    if row['intron_motif'] == 0 and row['unique_reads'] < 3:
        return False
    if row['intron_length'] > 50000 and row['unique_reads'] < 2:
        return False
    if row['intron_length'] > 100000 and row['unique_reads'] < 3:
        return False
    if row['intron_length'] > 200000 and row['unique_reads'] < 4:
        return False
    if row['intron_motif'] == 0 and row['max_overhang'] < 30:
        return False
    if row['intron_motif'] != 0 and row['max_overhang'] < 12:
        return False
    return True

# Apply the filter
df_filtered = df[df.apply(filter_junctions, axis=1)]

# Save the filtered file in the output folder
output_file = f"${output_folder}/${base_name}.SJ.out.filtered.tab"
df_filtered.to_csv(output_file, sep='\t', index=False, header=False)

print(f"Processed {file_name}: Original junctions: {len(df)}, Filtered junctions: {len(df_filtered)}, Removed: {len(df) - len(df_filtered)}")
END

done
