#!/bin/bash

# change directory
cd /research/labs/neurology/fryer/projects/sepsis/pig/Ecoli/PIPseq_snRNAseq/2024_Brain/usftp21.novogene.com/

# create file with list of R1 samples
awk '{print $2}' MD5.txt | grep _1.fq > R1Samples_snRNA.txt

sed -i 's/\_1.fq/_R1.fastq/g' R1Samples_snRNA.txt

# loops through list and print first line
touch sampleReadInfo_snRNA.txt
for sample in `cat R1Samples_snRNA.txt`; do
    #printf "${sample}\t"
    zcat ${sample} | head -1 >> sampleReadInfo_snRNA.txt	
done;

mv R1Samples_snRNA.txt  /research/labs/neurology/fryer/m239830/Ecoli_pigs/snRNAseq/scripts/R1Samples_snRNA.txt
mv sampleReadInfo_snRNA.txt /research/labs/neurology/fryer/m239830/Ecoli_pigs/snRNAseq/scripts/sampleReadInfo_snRNA.txt

cd /research/labs/neurology/fryer/m239830/Ecoli_pigs/snRNAseq/scripts/
paste -d "\t" R1Samples_snRNA.txt sampleReadInfo_snRNA.txt > sampleReadGroupInfo_snRNA.txt
rm R1Samples_snRNA.txt
rm sampleReadInfo_snRNA.txt

