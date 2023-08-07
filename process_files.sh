#!/bin/bash

# Run the get_pedigree.py script on the Parentage_8pops.csv file
python get_pedigree.py Parentage_8pops.csv

# Loop over all files ending in _output.txt
for file in *_output.txt
do
    # Extract part of the filename to use as the directory name
    dir_name=$(echo $file | cut -d "_" -f1)
    
    # Create a new directory (if it doesn't already exist)
    mkdir -p $dir_name
    
    # Move the file into the new directory
    mv $file $dir_name/

    # Extract the second line from the file and save it to sorted_bams
    awk 'NR==2{for(i=3; i<=NF; i++) print $i}' $dir_name/$file > $dir_name/sorted_bams
done
