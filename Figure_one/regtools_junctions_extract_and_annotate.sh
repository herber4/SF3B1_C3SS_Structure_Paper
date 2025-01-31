#!/bin/bash

regtools junctions extract -s 1 -a 8 -m 50 -M 500000 $B -o $N.junc

regtools junctions annotate K72.junc GRCh38.primary_assembly.genome.fa gencode.v39.annotation.gtf -o K72_junc_annotated.txt
