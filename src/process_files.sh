#!/bin/bash

# Run the get_pedigree.py 
python get_pedigree.py Parentage_8pops.csv

# Loop over all file
for file in *_output.txt
do
   
    dir_name=$(echo $file | cut -d "_" -f1)
    
    
    mkdir -p $dir_name

    mv $file $dir_name/

    awk 'NR==2{for(i=3; i<=NF; i++) print $i}' $dir_name/$file > $dir_name/sorted_bams
done
