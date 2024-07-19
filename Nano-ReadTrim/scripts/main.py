#!/usr/bin/env python3

import sys
import subprocess
from utils import find_sample_fastqs

query_loc = "/queries/"
outdir_loc = "/results/"
tool_choice = sys.argv[1]
gzip_output = sys.argv[2]
other_input_flags = sys.argv[3:]

# Convert list of flags into space-separated string
flat_flags = ""
for flag in other_input_flags:
    flat_flags = flat_flags + flag + " "
if len(flat_flags) > 0:
    flat_flags = flat_flags[:-1]

# Find samples (and PE FASTQ locs) within the fastq_dir
sampleid_fastq_dict = find_sample_fastqs(query_loc)

for base_id in sampleid_fastq_dict:
    R1, R2 = sampleid_fastq_dict[base_id]

    if tool_choice == "fastp":
        command = "fastp" + \
                    " -i " + R1 + \
                    " -I " + R2
        if gzip_output == "TRUE":
            command += " -o " + outdir_loc + base_id + "-trim_R1.fastq.gz" + \
                        " -O " + outdir_loc + base_id + "-trim_R2.fastq.gz"
        elif gzip_output == "FALSE":
            command += " -o " + outdir_loc + base_id + "-trim_R1.fastq" + \
                        " -O " + outdir_loc + base_id + "-trim_R2.fastq"
        else:
            raise Exception("Unallowed gzip_out (TRUE or FALSE):", gzip_output)
        
    elif tool_choice == "trimmomatic":
        if gzip_output == "TRUE":
            forward_paired = outdir_loc + base_id + "-trim_R1.fastq.gz"
            reverse_paired = outdir_loc + base_id + "-trim_R2.fastq.gz"
            forward_unpaired = "/dev/null/R1.fastq.gz"
            reverse_unpaired = "/dev/null/R2.fastq.gz"
        elif gzip_output == "FALSE":
            forward_paired = outdir_loc + base_id + "-trim_R1.fastq"
            reverse_paired = outdir_loc + base_id + "-trim_R2.fastq"
            forward_unpaired = "/dev/null/R1.fastq"
            reverse_unpaired = "/dev/null/R2.fastq"
        else:
            raise Exception("Unallowed gzip_out (TRUE or FALSE):", gzip_output)

        command = "TrimmomaticPE" + \
                   " " + R1 + " " + R2 + \
                   " " + forward_paired + " " + forward_unpaired + \
                   " " + reverse_paired + " " + reverse_unpaired

    command += " " + flat_flags
    process = subprocess.run(command.split(), stdout=subprocess.PIPE)
    if process.returncode != 0:
        raise Exception("Error running read trimming", command, "Error:", process.returncode)
