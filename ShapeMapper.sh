#!/bin/bash

ml shapemapper/2.2.0

for f in *.fasta; do
    N=$(basename ${f} .fasta);
    shapemapper --nproc 16 --name $N --out $N --target $f --modified --R1 ../e1/* --R2 ../e2/* --untreated --R1 ../c1/* --R2 ../c2/* --overwrite ;
done
