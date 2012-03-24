#!/bin/bash

python test.py test_gen.sam
samtools view -bS test_gen.sam > test_gen.bam
~/bamtools/bin/bamtools sort -in test_gen.bam -out test_gen_sorted.bam -byname
../bobby -o combined.bam test_gen_sorted.bam 
samtools view combined.bam 

#rm combined.bam
#rm test_gen.bam
#rm test_gen.sam
#rm test_gen_sorted.bam
