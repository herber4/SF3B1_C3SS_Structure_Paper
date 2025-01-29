

module load hisat2/2.1.0
module load samtools/1.10

for i in `ls -1 *_clean_R1.fastq.gz | sed 's/_clean_R1.fastq.gz//'`
do
hisat2 -p 24 --qc-filter -x /data2/lackey_lab/austin/dbs/spike_mas/hisat_index -1 ${i}_clean_R1.fastq.gz -2 ${i}_clean_R2.fastq.gz -S ${i}.sam
done


for B in *.sam; do
	N=$(basename $B .sam) ;
	samtools view --threads 8 -bS $B | samtools sort -O BAM -o $N.bam ;
done

for B in *.bam; do
  samtools index $B ;
done
