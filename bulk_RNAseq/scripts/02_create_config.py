#!/usr/bin/python3

# create a new output file
outfile = open('config.json', 'w')

# get all, kidney, brain and blood sample names
allSamples = list()
numSamples = 0

with open('sampleReadGroupInfo.txt', 'r') as infile:
    for line in infile:
        numSamples += 1

        line = line.replace(".", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # E1_BR.FCHVC2VDRXY_L1_R1_ITAAGTGGT-CTTAAGCC.fastq.gz pigID_tissue.sequencer_lane_read_X-X.fastq.gz

        # create a shorter sample name
        stemName = sampleAttributes[0] + '_' + sampleAttributes[1]
        allSamples.append(stemName)

# create header and write to outfile
header = '''{{
    "Commment_Input_Output_Directories": "This section specifies the input and output directories for scripts",
    "rawReads" : "/research/labs/neurology/fryer/projects/sepsis/pig/Ecoli/bulkRNA/",
    "rawQC" : "../rawQC/",
    "trimmedReads" : "../trimmedReads/",
    "trimmedQC" : "../trimmedQC/",
    "starAligned" : "../starAligned/",
    "bamstats" : "../bamstats/",
    "kallisto" : "../kallisto/",


    "Comment_Reference" : "This section specifies the location of the Sus scrofa, Ensembl reference genome",
    "Sscrofa.Ymasked.fa" : "/research/labs/neurology/fryer/projects/references/pig/ensembl_v7/Sus_scrofa.Sscrofa11.1.dna.toplevel.Ymask_2024",
    "Scrofa.cdna.Ymasked.fa" : "/research/labs/neurology/fryer/projects/references/pig/ensembl_v7/Sus_scrofa.Sscrofa11.1.cdna.all.Ymask",
    "Sscrofa.fa" : "/research/labs/neurology/fryer/projects/references/pig/ensembl_v7/Sscrofa11.1.dna.toplevel",
    "Sscrofa.gtf" : "/research/labs/neurology/fryer/projects/references/pig/ensembl_v7/Sus_scrofa.Sscrofa11.1.107",
    "star_ref_index" : "/research/labs/neurology/fryer/projects/references/pig/ensembl_v7/Sus_scrofa.Sscrofa11.1.dna.toplevel_star_Ymask",
    "kallisto_ref_index" : "/research/labs/neurology/fryer/projects/references/pig/ensembl_v7/Sus_scrofa.Sscrofa11.1.cdna.all.Ymask.kallisto",

    "Comment_Sample_Info": "The following section lists the samples that are to be analyzed",
    "sample_names": {0},
'''
outfile.write(header.format(allSamples))

# config formatting
counter = 0
with open('sampleReadGroupInfo.txt', 'r') as infile:
    for line in infile:
        counter += 1
        # store sample name and info from the fastq file
        split = line.split()
        base = split[0]
        base = base.replace(".fastq.gz", "")
        sampleName1 = base
        sampleName2 = sampleName1.replace("R1","R2")
        sampleName3 = base.replace("L1","L2")
        sampleName4 = sampleName2.replace("L1","L2")
        base = base.replace("_R1_001", "")
        sampleInfo = split[1]

        # make naming consistent, we will rename using only underscores (no hyphens)
        line = line.replace(".", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')  # uniqueNum_pigID_tissue._XX_XX_sequencer_lane_read_001.fastq.gz

        # create a shorter sample name
        stemName = sampleAttributes[0] + '_' + sampleAttributes[1] 
        shortName1 = stemName + '_R1'
        shortName2 = stemName + '_R2'

        # break down fastq file info
        # @A00127:312:HVNLJDSXY:2:1101:2211:1000
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
        sampleInfo = sampleInfo.split(':')
        instrument = sampleInfo[0]
        runNumber = sampleInfo[1]
        flowcellID = sampleInfo[2]

        lane = sampleInfo[3]
        ID = stemName  # ID tag identifies which read group each read belongs to, so each read group's ID must be unique
        SM = stemName  # Sample
        PU = flowcellID + "." + lane  # Platform Unit
        LB = stemName

        out = '''
    "{0}":{{
        "fq_path": "/research/labs/neurology/fryer/projects/sepsis/pig/Ecoli/bulkRNA/bulkRNA_mereged_lanes/",
        "fq1": "{1}",
        "fq2": "{2}",
        "fq3": "{3}",
        "fq4": "{4}",
        "fq_R1": "{5}",
        "fq_R2": "{6}",
        "ID": "{7}",
        "SM": "{7}",
        "PU": "{8}",
        "LB": "{9}",
        "PL": "Illumina"
        '''
        outfile.write(out.format(stemName, sampleName1, sampleName2, sampleName3, sampleName4, shortName1, shortName2, stemName, PU, LB))
        if (counter == numSamples):
            outfile.write("}\n}")
        else:
            outfile.write("},\n")
outfile.close()
