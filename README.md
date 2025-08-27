# Ecoli_pigs
Bulk and single nucleus RNAseq of Sus scrofa (pigs) that received either saline or Escherichia coli  (E. coli).

The goal of this experiment is to identify differentially expressed genes (DEGs) between experimental groups.  Pigs were injected with saline (control) or Escherichia coli  (E. coli) to model sepsis.  Brain cotrex samples were collected and sent for bulk RNA and single nucleus RNA sequencing.

Explore gene expresison within the single nucleus data in our published [shiny app](https://fryerlab.shinyapps.io/Ecoli_snRNAseq/)


## Set up conda environment
This workflow uses conda. For information on how to install conda [here](https://docs.conda.io/projects/conda/en/latest/user-guide/index.html)

To create the environment:
```
conda env create -n Ecoli_pigs --file Ecoli_pigs.yml # for bulk RNAseq
conda env create -n cellbender_py38 --file cellbender_py38.yml # for single nucleus RNAseq

# To activate this environment, use
#
#     $ conda activate Ecoli_pigs
#
# To deactivate an active environment, use
#
#     $ conda deactivate Ecoli_pigs

```
## Bulk RNAseq differential expression
We have put together a workflow for inferring differential expression between E. coli and saline (control) female pigs using read aligner STAR. The tools used in this workflow are publicly available and we ask that if you use this workflow to cite the tools used. 

### 1. Download fastq files and pig reference genome. Create a reference index. 
The raw fastq files may be obtained from SRA PRJNA1172687 - bulk; PRJNA1175698 - snRNAseq. Bulk RNAseq samples were sequenced to ~50 million (M) 2 × 100 bp paired-end reads across two lanes. Single nucleus RNAseq aimed to capture 15,000 nuclei per sample. Information on how to download from SRA may be found [here](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/). 

Download the Sus scrofa (pig) reference genome and gene annotation from Ensembl. The version used in this workflow is v7. 
```
cd reference
wget http://ftp.ensembl.org/pub/release-107/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz
wget http://ftp.ensembl.org/pub/release-107/fasta/sus_scrofa/cdna/Sus_scrofa.Sscrofa11.1.cdna.all.fa.gz
wget http://ftp.ensembl.org/pub/release-107/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.103.gtf.gz 
```

All pigs used in this study are genetically XX female. To avoid mis-mapping of homologous X and Y-linked genes, samples were aligned to a reference with the Y chromosome hard masked. See [Olney et al. 2020](https://bsd.biomedcentral.com/articles/10.1186/s13293-020-00312-9) for more details about this approach. 

Build the reference genome index:
```
python hardmaskY.py

STAR --runThreadN 12 --runMode genomeGenerate --genomeDir Sus_scrofa.Sscrofa11.1.dna.toplevel_star_Ymask --genomeFastaFiles Sus_scrofa.Sscrofa11.1.dna.toplevel.Ymask.fa --sjdbGTFfile Sus_scrofa.Sscrofa11.1.107.gtf
```

### 2. Align reads and generate quantification estimates.
First move to the scripts snakemake folder.
```
cd bulkRNA/scripts/snakemake/
```
Now run the snakefile. To run multiple samples in parallel use -j and the number of jobs to run.
You may need to adjust pathways within the config.json file to be to the location of the files on your system. 
```
snakemake -s Snakefile
```

### 3. preform bulk differntial expression
First, you will need to create sub-directories within the results folders. This is where the results from running differential expression will be placed. 
```
cd results/
mkidr CPM  DEGs  JSD  MDS  boxplot  density  gprofiler  metascape  volcano  voom
```

Move to the scripts R folder.
```
cd ../scripts/R/
R 04a_stats.R # Statistical analysis of clinical and pathological measurements
R 04_differential_expression_4Ecoli_vs_4Saline.Rmd # gene-level differential expression 
```


### 4. single nuclues RNAseq data alignment, QC, and differential expression. 
Frozen brain samples from the same pigs used for bulk RNAseq analysis (n = 4 pigs per group) were used for  single nucleus. Alignment of the single nucleus 10x Genomics data to the Ensembl reference genome Sscrofa11.1 v107 was carried out using the Cellranger. Cellbender v0.3.2 was used to remove technical noise from the single nucleus data that may result from enzymatic processes and produce library fragments leading to contamination.

```
cd /snRNAseq/scripts/

conda activate cellbender_py38

# create config file for snakemake pipeline to align the reads and remove ambient RNA
sh 00_rename_fastq_files.hs
sh 01_get_read_info.sh
py 02_create_config.py
sh 03_run_Snakefile.sh # must update to your HPC system
```

The counts adjusted single nucleus count data from Cellbender was then imported into R v4.3.0 using the Seurat v5.0.2 package for quality filtering, dimensionality reduction, clustering, QC, differential expression, and cellchat analysis. R scripts 04 through 21. 

## Contacts

| Contact | Email |
| --- | --- |
| Kimberly Olney, PhD | kolney@tgen.org |
| John Fryer, PhD | jfryer@tgen.org |
