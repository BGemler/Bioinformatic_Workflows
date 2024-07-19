#!/usr/bin/env python3

import sys
import os
import subprocess
from utils import find_sample_fastqs, search_for_threads_in_flags
from multiprocessing import Process, Queue, Pool, Manager

tool_choice = sys.argv[1]
gzip_output = sys.argv[2]
fastq_dir = "/queries/"
out_dir = "/results/"
ref_dir = "/genome_dir/" + sys.argv[3]

other_input_flags = sys.argv[4:]

# Convert list of flags into space-separated string
flat_flags = ""
for flag in other_input_flags:
    flat_flags = flat_flags + flag + " "
if len(flat_flags) > 0:
    flat_flags = flat_flags[:-1]

# Find samples (and PE FASTQ locs) within the fastq_dir
sampleid_fastq_dict = find_sample_fastqs(fastq_dir)

for base_id in sampleid_fastq_dict:
    R1, R2 = sampleid_fastq_dict[base_id]

    if tool_choice == "bowtie2":
        command = tool_choice + \
                    " -i " + R1 + \
                    " -I " + R2 + \
                    " -x " + ref_dir
        if gzip_output == "TRUE":
            command += " --un-conc-gz " + out_dir + base_id + "-hostr_R%%.fastq.gz"
        elif gzip_output == "FALSE":
            command += " --un-conc " + out_dir + base_id + "-hostr_R%%.fastq"
        else:
            raise Exception("Unallowed gzip_output (required TRUE or FALSE):", gzip_output)

    elif tool_choice == "kraken2":
        command = tool_choice + \
                    " --paired " + R1 + " " + R2 + \
                    " --db " + ref_dir + \
                    " --unclassified-out " + out_dir + base_id + "-hostr_R#.fastq"
    
    else:
        raise Exception("Unallowed tool_choice (required kraken2 or bowtie2):", tool_choice)

    # add flat flags
    command += " " + flat_flags

    # execute host removal
    process = subprocess.run(command.split(), stdout=subprocess.PIPE)
    if process.returncode != 0:
        raise Exception("Error running host removal", command, "Error:", process.returncode)

    # If kraken2 is run, we want to correct the output to be _R1, NOT _1
    if tool_choice == "kraken2":
        command = "mv " + out_dir + base_id + "-hostr_R_1.fastq" + " " + out_dir + base_id + "-hostr_R1.fastq"
        process = subprocess.run(command.split(), stdout=subprocess.PIPE)
        if process.returncode != 0:
            raise Exception("Error moving k2 output file", process.returncode)     

        command = "mv " + out_dir + base_id + "-hostr_R_2.fastq" + " " + out_dir + base_id + "-hostr_R2.fastq"
        process = subprocess.run(command.split(), stdout=subprocess.PIPE)
        if process.returncode != 0:
            raise Exception("Error moving k2 output file", process.returncode)         


def run_command(command):
	os.system(command)
	return

# If kraken2 is tool and Gzip = TRUE, then gzip each .fastq file
# Note - do this at the end for parallelization advantages 
if tool_choice == "kraken2" and gzip_output == "TRUE":
    # This is very slow since it's single threaded!
    # We want to multi-thread this process. 
    num_threads = search_for_threads_in_flags(other_input_flags)
    print("Gzipping Kraken2 output using # threads:", num_threads)

    # Initiate multiprocessing
    m = Manager()
    q = m.Queue()
    p = Pool(int(num_threads))

    commands_to_run = []
    for base_id in sampleid_fastq_dict:
        file_locs = [out_dir + base_id + "-hostr_R1.fastq", \
                    out_dir + base_id + "-hostr_R2.fastq"]
        for file_loc in file_locs:
            command = "gzip " + file_loc
            commands_to_run.append(p.apply_async(run_command, (command,)))

    [r.get() for r in commands_to_run]
