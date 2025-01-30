#!/bin/bash

ml viennarna/2.4.18

for f in *.fasta; do
    # Extract the base name without extension
    N=$(basename "$f" .fasta)
    
    # Extract the gene name (before "_md_")
    gene_name=$(echo "$N" | sed 's/_md_.*//')
    
    # Construct the shape file path
    shape_file="../shape/${gene_name}.shape"
    
    # Extract the threshold number after "_md_"
    threshold=$(echo "$N" | sed 's/.*_md_//')
    
    # Run RNAfold with the dynamically constructed shape file
    RNAfold -p --bppmThreshold=0 --maxBPspan $threshold --shape ../shape/"$shape_file" --noPS "$f" > "$N.out"
done


cp *.ps base_pairs/
cd base_pairs

./extract_base_probs.sh
