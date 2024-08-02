#!/usr/bin/env python3

import sys
import os
import subprocess

# Set temporary matplotlib directory to aviod warning message
os.environ['MPLCONFIGDIR'] = "/tmp/"
import matplotlib.pyplot as plt

from utils import find_sample_fastqs

spades_exe = "/SPAdes-4.0.0-Linux/bin/spades.py"
quast_exe = "python3 /quast-5.2.0/quast.py"

query_loc = "/queries/"
outdir_loc = "/results/"

refgenome_loc = sys.argv[1].split("/")[-1]
refgene_loc = sys.argv[2].split("/")[-1]

other_input_flags = sys.argv[3:]

# Make plotting output loc
plotting_out_loc = outdir_loc + "contig-plots/"
os.makedirs(plotting_out_loc, exist_ok=True)

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

    base_id_result_dir = outdir_loc + base_id + "/"
    #os.makedirs(base_id_result_dir, exist_ok=True)

    command = spades_exe + \
                " -1 " + R1 + \
                " -2 " + R2 + \
                " -o " + base_id_result_dir
    command += " " + flat_flags
    print(command)
    os.system(command)
    """
    process = subprocess.run(command.split(), stdout=subprocess.PIPE)
    if process.returncode != 0:
        raise Exception("Error running SPAdes", command, "Error:", process.returncode)
    """

    # get contigs, coverage
    contig_loc = base_id_result_dir + "contigs.fasta"
    all_coverages = []
    with open(contig_loc, "r") as f:
        for row in f:
            if row.startswith(">"):
                coverage = float(row.split("_")[-1])
                all_coverages.append(coverage)
    f.close()

    # Make histogram of the coverage values
    plt.hist(all_coverages, bins = 20)
    plt.title(base_id + " # of Contigs:" + str(len(all_coverages)))
    plt.xlabel("Coverage Value")
    plt.ylabel("# of Contigs per Bin")
    plt.savefig(plotting_out_loc + base_id + "-contig_coverages.png", bbox_inches = 'tight')
    plt.close()

    # Run QUAST to evaluate assembly quality 
    quast_out_loc = base_id_result_dir + "quast_output/"
    command = quast_exe + \
                " " + contig_loc + \
                " -1 " + R1 + \
                " -2 " + R2 + \
                " -o " + quast_out_loc
    if refgenome_loc != "NA":
        command += " -r " + "/refgenomes/" + refgenome_loc
    if refgene_loc != "NA":
        command += " -g " + "/refgenomes/" + refgene_loc
    
    os.system(command)
