#load mod
module load bbmap/38.73
module load java_jdk/11.0.2

for i in `ls -1 *_R1.fastq.gz | sed 's/_R1.fastq.gz//'`
do
bbduk.sh -Xmx64g in1=${i}_R1.fastq.gz in2=${i}_R2.fastq.gz out1=${i}_clean_R1.fastq.gz out2=${i}_clean_
R2.fastq.gz ref=/data2/databases/rrna_silva/ribokmers.fa,adapters.fa tbo tpe ktrim=r k=31 refstats=$i.t
xt hdist=1 qtrim=rl trimq=15 minlength=75 maxns=1 threads=auto
done
