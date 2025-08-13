# Ecoli_pigs
Bulk RNAseq of Sus scrofa (pigs) that received either saline or Escherichia coli  (E. coli).

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
We have put together a workflow for inferring differential expression between saline (control) and LPS female pigs using two read aligners STAR and Kallisto. These tools are publicly available and we ask that if you use this workflow to cite the tools used listed in the table below. 

### 1. Download fastq files and pig reference genome. Create a reference index. 
The raw fastq files may be obtained from SRA PRJNA723823. Samples were sequenced to ~50 million (M) 2 × 100 bp paired-end reads across two lanes. Information on how to download from SRA may be found [here](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/). 

Download the Sus scrofa (pig) reference genome, transcriptome, and gene annotation from Ensembl. The version used in this workflow is v7. 
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

Build the reference transcriptome index:

```
pyhon hardmaskY_transcriptome.py

kallisto index -i Sus_scrofa.Sscrofa11.1.cdna.all.Ymask.kallisto.fa Sus_scrofa.Sscrofa11.1.cdna.all.Ymask.fa
```

### 2. Align reads and generate quantification estimates.
First move to the scripts snakemake folder.
```
cd bulkRNA/scripts/snakemake/
```
Now run the snakefile. To run multiple samples in parallel use -j and the number of jobs to run.
```
snakemake -s Snakefile
```

### 3. preform differntial expression
First, you will need to create sub-directories within the star and kallisto folders found in the results directory. This is where the results from running differential expression will be placed. 
```
cd results/star/
mkidr CPM  DEGs  JSD  MDS  boxplot  density  gprofiler  metascape  volcano  voom
cd ../kallisto/
mkdir CPM  DEGs  JSD  MDS  boxplot  density  gprofiler  metascape  volcano  voom TPM
```
Move to the scripts R folder.
```
cd ../../scripts/R/
R gene_DE.R # gene-level differential expression 
R pseudoalign_isoform_DE.Rmd # isoform-level differential expression 
R LPS_mouse_processing.Rmd # processing of the mouse Kang et al. 2017 data
```

Additional R scripts are for comparing between tissues within the pig data, and between mouse and pig brain. Scripts for making the manuscript figures are located in the folder called manuscript_figures.

## Contacts

| Contact | Email |
| --- | --- |
| Kimberly Olney, PhD | olney.kimberly@mayo.edu |
| John Fryer, PhD | fryer.john@mayo.edu |
