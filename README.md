# Ecoli_pigs
Bulk and single nucleus RNAseq of Sus scrofa (pigs) that received either saline or Escherichia coli  (E. coli); n = 4 per group. 

The goal of this experiment is to identify differentially expressed genes (DEGs) between experimental groups.  Pigs were injected with saline (control) or Escherichia coli  (E. coli) to model sepsis.  Brain cotrex samples were collected and sent for bulk RNA and single nucleus RNA sequencing.

Explore gene expresison in our published shiny apps: [single nucleus](https://fryerlab.shinyapps.io/Ecoli_snRNAseq/) and [bulk RNAseq](https://fryerlab.shinyapps.io/Ecoli_snRNAseq/)

| Sample     | Group   | 
| -----------|:-------:|
| E1           | E. coli  | 
| E4           | E. coli  | 
| E6           | E. coli  | 
| E8           | E. coli  | 
| S1           | Saline   | 
| S2           | Saline   | 
| S3           | Saline   | 
| S6           | Saline   | 


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
The raw fastq files may be obtained from SRA PRJNA1172687 - bulk; PRJNA1175698 - snRNAseq. Please note that the SRA bulk project contains kidney tissue as well, which was not used in this study. Bulk RNAseq samples were sequenced to ~50 million (M) 2 × 100 bp paired-end reads across two lanes. Single nucleus RNAseq aimed to capture 15,000 nuclei per sample. Information on how to download from SRA may be found [here](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/). 

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
cd bulk_RNAseq/scripts/

# The config file that is needed to run snakemake was created using scripts 1-2. 
# obtain sample information from the fastq files
sh 01_get_readgroup_info.sh 
# output is sampleReadGroupInfo.txt

# builds the config file from the output of 01_get_readgroup_info
python 02_create_config.py
# output is config.json
```

Now run the snakefile. To run multiple samples in parallel use -j and the number of jobs to run.
You may need to adjust pathways within the config.json file to be to the location of the files on your system. 
```
snakemake -s Snakefile
```
The output is STAR aligned bam files and counts information

### 3. preform bulk differntial expression
This section contains R mark down files for comparing between groups for clinical and pathological measurements and gene-level differential expression, followed by GO enrichment analysis. 
```
R 04a_stats.Rmd # Statistical analysis of clinical and pathological measurements
R 04_differential_expression_4Ecoli_vs_4Saline.Rmd # gene-level differential expression
R 05_gprofiler_GO_terms.Rmd # Gene ontology form gProfiler
```
Additionally, metascape was used to obtain GO enrichment terms, for this the web based version of metascape was used. The results were saved in results/metascape and used to generate the GO enrichment plots.  


## Single-Nucleus RNA Sequencing (snRNA-seq) Workflow
To resolve cell-type-specific transcriptional responses, we performed single-nucleus RNA sequencing (snRNA-seq) on frozen cortical tissue from the same experimental cohort (n=4 per group). Data Processing and Quality ControlRaw sequencing reads were aligned to the Ensembl Sscrofa11.1 (v107) reference genome using Cell Ranger (10x Genomics). To mitigate technical noise and ambient RNA contamination—often exacerbated by enzymatic dissociation in porcine brain tissue—we employed CellBender (v0.3.2) for background removal.\
Downstream analysis was conducted in R using the Seurat framework. Potential doublets were identified and excluded using DoubletFinder. Major brain cell populations were identified and annotated based on the expression of canonical lineage markers Differential Expression and Pathway AnalysisWe utilized a pseudobulk approach to identify differential expression between E. coli and saline control groups within each distinct cell type. Functional enrichment was characterized using Metascape to identify significantly altered Gene Ontology (GO) terms.\
Repository Structure and Implementation:\
The analysis pipeline is contained within the /snRNAseq/scripts/ directory.
Input Data: Raw FASTQ files should be located in snRNAseq/fastq/
Reference Genome: Scripts point to ../../projects/references/pig/ensembl_v7.
**Note**: Update the paths within the scripts to match your local environment configuration before execution.

```
cd /snRNAseq/scripts/
conda activate cellbender_py38 # this contains the cellbender tool
```

### create config file for snakemake pipeline to align the reads and remove ambient RNA
```
sh 00_rename_fastq_files.sh
sh 01_get_read_info.sh
py 02_create_config.py
sh 03_run_Snakefile.sh # must update to your HPC job submission system
```

The counts adjusted single nucleus count data from Cellbender was then imported into R v4.3.0 using the Seurat v5.0.2 package for quality filtering, dimensionality reduction, clustering, QC, differential expression, and cellchat analysis. R scripts 04 through 21. 

### Scripts used to make the figures for the manuscript are located in manuscript_figure_scripts
```
```

## Contacts

| Contact | Email |
| --- | --- |
| Kimberly Olney, PhD | kolney@tgen.org |
| John Fryer, PhD | jfryer@tgen.org |
