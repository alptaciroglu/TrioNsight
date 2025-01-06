#!/bin/bash

# Directory containing the dat and pdb files
INPUT_DIR="../models_aligned/"

# Directory for output files
OUTPUT_DIR="../strucclus_res/"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all dat files in the input directory
for dat_file in "$INPUT_DIR"/*.dat; do
    # Extract the base name without extension
    base_name=$(basename "$dat_file" .dat)
    
    # Corresponding pdb file
    pdb_file="$INPUT_DIR/$base_name.pdb"
    
    # Check if pdb file exists
    if [ ! -f "$pdb_file" ]; then
        echo "Error: Corresponding PDB file not found for $dat_file"
        continue
    fi
    
    # Get the number of rows in the dat file
    num_rows=$(wc -l < "$dat_file")
    
    # Determine the upper limit for n (maximum of num_rows or 8, whichever is lower)
    upper_limit=$((num_rows < 10 ? num_rows : 10))
    
    # Inner loop for n values
    for ((n=1; n<=upper_limit; n++)); do
        # Generate output file name
        output_file="$OUTPUT_DIR/${base_name}_n${n}.txt"
        
        # Run the command and redirect output to the file
        ./strucclus -n $n -t average "$dat_file" "$pdb_file" > "$output_file"
        
        echo "Processed $dat_file with n=$n, output written to $output_file"
    done
done