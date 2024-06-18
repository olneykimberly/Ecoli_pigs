#!/bin/bash

# change directory
cd /research/labs/neurology/fryer/projects/sepsis/pig/Ecoli/bulkRNA/

# create file with list of R1 samples
ls -1 | grep L1_R1_ > L1R1Samples.txt

# change directory 

# loops through list 
touch sampleReadInfo.txt
for sample in `cat L1R1Samples.txt`; do
    zcat ${sample} | head -1 >> sampleReadInfo.txt
done;

# mv the files 
mv L1R1Samples.txt  /research/labs/neurology/fryer/m239830/Ecoli_pigs/bulk_RNAseq/scripts/L1R1Samples.txt
mv sampleReadInfo.txt /research/labs/neurology/fryer/m239830/Ecoli_pigs/bulk_RNAseq/scripts/sampleReadInfo.txt

cd /research/labs/neurology/fryer/m239830/Ecoli_pigs/bulk_RNAseq/scripts/
paste -d "\t" L1R1Samples.txt sampleReadInfo.txt > sampleReadGroupInfo.txt
rm L1R1Samples.txt
rm sampleReadInfo.txt
