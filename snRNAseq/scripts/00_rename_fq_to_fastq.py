import os
import shutil

def rename_files_in_directories(directories):
    for directory in directories:
        directory_path = os.path.join("/research/labs/neurology/fryer/projects/sepsis/pig/Ecoli/PIPseq_snRNAseq/2024_Brain/usftp21.novogene.com/01.RawData/", directory)
        for root, dirs, files in os.walk(directory_path):
            for file in files:
                if file.endswith("_1.fq.gz"):
                    old_path = os.path.join(root, file)
                    new_path = os.path.join(root, file.replace("_1.fq.gz", "_R1.fastq.gz"))
                    shutil.move(old_path, new_path)
                    print(f"Renamed: {old_path} -> {new_path}")
                elif file.endswith("_2.fq.gz"):
                    old_path = os.path.join(root, file)
                    new_path = os.path.join(root, file.replace("_2.fq.gz", "_R2.fastq.gz"))
                    shutil.move(old_path, new_path)
                    print(f"Renamed: {old_path} -> {new_path}")

# List of directories to process
directories_to_process = ["E2A", "E2B", "Undetermined"]

rename_files_in_directories(directories_to_process)

