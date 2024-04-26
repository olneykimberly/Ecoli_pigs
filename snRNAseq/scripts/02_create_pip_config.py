#!/usr/bin/python3

# create a new output file
outfile = open('config.json', 'w')

allSamples = list()
sex = list()
read = ["R1", "R2"]
numSamples = 0

with open('sampleReadGroupInfo_snRNA.txt', 'r') as infile:
    for line in infile:
        numSamples += 1

       # line = line.replace(".", "_")
        line = line.replace("/", "_")
        split = line.split()
        sampleAttributes = split[0].split('_') # E2A/E2A_CKDL240011999-1A_22JTHGLT3_L1_R1.fastq.gz
        # create a shorter sample name
        stemName = sampleAttributes[1]
        allSamples.append(stemName)

# create header and write to outfile
header = '''{{
    "Commment_Input_Output_Directories": "This section specifies the input and output directories for scripts",
    "pipseeker_dir" : "../pipseeker_output/", 
    "cellbender_dir" : "../cellbender/",
    "results" : "../results/",
    "rObjects" : "../rObjects/",
    "fastq_path" : "/research/labs/neurology/fryer/projects/sepsis/pig/Ecoli/PIPseq_snRNAseq/2024_Brain/usftp21.novogene.com/01.RawData/",

    "Comment_Reference" : "This section specifies the location of the Sus scrofa, Ensembl reference genome",
    "Sus.scrofa" : "/research/labs/neurology/fryer/projects/references/pig/ensembl_v7/Sus_scrofa.Sscrofa11.1.dna.toplevel_star_Ymask_filtered_gtf/",

    "Comment_Sample_Info": "The following section lists the samples that are to be analyzed",
    "sample_names": {0},
    "read": {1},
'''
outfile.write(header.format(allSamples, read))

# config formatting
counter = 0
with open('sampleReadGroupInfo_snRNA.txt', 'r') as infile:
    for line in infile:
        counter += 1
        # store sample name and info from the fastq file
        split = line.split()
        base = split[0]
        base = base.replace("_R1.fastq.gz", "")
        sampleName1 = base
        sampleInfo = split[1]

        # make naming consistent, we will rename using only underscores (no hyphens)
        line = line.replace("/", "_")
        split = line.split()
        sampleAttributes = split[0].split('_')
        # uniqueNum-number_sequencer_lane_read.fastq.gz

        # create a shorter sample name
        stemName = sampleAttributes[1]
        stemID = sampleAttributes[1]
        fullName = sampleAttributes[0] + '_' + sampleAttributes[1] + '_' + sampleAttributes[2] + '_' + sampleAttributes[3] 
        shortName1 = stemName + '_R1'
        shortName2 = stemName + '_R2'

        # break down fastq file info
        # @A00127:312:HVNLJDSXY:2:1101:2211:1000
        # @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>
        sampleInfo = sampleInfo.split(':')
        instrument = sampleInfo[0]
        runNumber = sampleInfo[1]
        flowcellID = sampleInfo[2]

        lane = sampleInfo[6]
        ID = stemID  # ID tag identifies which read group each read belongs to, so each read group's ID must be unique
        SM = fullName  # Sample
        PU = flowcellID + "." + lane  # Platform Unit
        LB = stemName

        out = '''
    "{0}":{{
        "fq_path": "/research/labs/neurology/fryer/projects/sepsis/pig/Ecoli/PIPseq_snRNAseq/2024_Brain/usftp21.novogene.com/",
        "fq": "{1}",
        "shortName1": "{2}",
        "shortName2": "{3}",
        "ID": "{4}",
        "SM": "{4}",
        "PU": "{5}",
        "LB": "{6}",
        "PL": "Illumina"
        '''
        outfile.write(out.format(stemName, sampleName1, shortName1, shortName2, stemID, PU, LB))
        if (counter == numSamples):
            outfile.write("}\n}")
        else:
            outfile.write("},\n")
outfile.close()
