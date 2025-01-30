#!/bin/bash
for f in *.fasta; do
	N=$(basename ${f} .fasta);
	python /ScanFold/ScanFold.py $f -s 20 -w 200 --react ../shape/${N}.shape --out_name $N ;
done
