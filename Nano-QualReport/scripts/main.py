#!/usr/bin/env python3

import sys
import subprocess
import os
from multiprocessing import Process, Queue, Pool, Manager

def run_command(command):
	os.system(command)
	return

# options/paths
number_threads = sys.argv[1]
query_loc = "/queries/"
outdir_loc = "/results/"

# Initiate multiprocessing
m = Manager()
q = m.Queue()
p = Pool(int(number_threads))

# Add a FastQC subdirectory 
fastqc_dir = outdir_loc + "fastqc_rpts/"
os.makedirs(fastqc_dir, exist_ok=True)

# run FastQC over all FASTQs in directory
commands_to_run = []
for filename in os.listdir(query_loc):
    command = "fastqc" + \
                " " + query_loc + filename + \
                " --outdir=" + fastqc_dir + \
                " -t 1 --quiet"
    commands_to_run.append(p.apply_async(run_command, (command,)))

# Execute all FastQC runs that are Q'd
[r.get() for r in commands_to_run]

# Run MultiQC over all FASTQC reports in directory
command = "multiqc" + \
            " " + fastqc_dir + \
            " -o " + outdir_loc + \
            " -f"
process = subprocess.run(command.split(), stdout=subprocess.PIPE)
if process.returncode != 0:
	raise Exception("Error running MultiQC", process.returncode)
