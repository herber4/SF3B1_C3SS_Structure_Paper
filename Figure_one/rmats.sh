
conda activate rmats
rmats.py --b1 control_bam_files.txt --b2 sf3b1_MT_bam_files.txt --gtf gencode.v39.annotation.gtf --task both -t paired --nthread 32 --novelSS --variable-read-length --readLength 150 --od N6_WT_controls_rMATS --tmp rMATS_p1_tmp
