#!/usr/bin/env python3

import sys
import os
import subprocess
from utils import find_sample_fastqs, create_feature_table, \
                    search_for_threads_in_flags, run_command, cal_alpha_diversity
from multiprocessing import Process, Queue, Pool, Manager

restrict_taxids_to_mb = sys.argv[1]

taxtree_loc = "/taxdmp/"
fastq_dir = "/queries/"
out_dir = "/results/"
nuc_ref_dir = "/nuc_index/"
prot_ref_dir = "/prot_index/"

bracken_exe = "/src/Bracken-2.9/bracken"

other_input_flags = sys.argv[2:]

# Convert list of flags into space-separated string
flat_flags = ""
for flag in other_input_flags:
    flat_flags = flat_flags + flag + " "
if len(flat_flags) > 0:
    flat_flags = flat_flags[:-1]

# Bracken is slow since it's single threaded!
# We want to multi-thread this process. 
num_threads = search_for_threads_in_flags(other_input_flags)
print("Running Bracken using # threads:", num_threads)
bracken_commands_to_run = []

m = Manager()
q = m.Queue()
p = Pool(int(num_threads))

# Find samples (and PE FASTQ locs) within the fastq_dir
sampleid_fastq_dict = find_sample_fastqs(fastq_dir)

for base_id in sampleid_fastq_dict:
    R1, R2 = sampleid_fastq_dict[base_id]

    # create sample-specific output folder
    base_id_out_fldr = out_dir + base_id + "/"
    os.makedirs(base_id_out_fldr, exist_ok=True)

    # Set up nucleotide command
    command = "kraken2" + \
                " --paired " + R1 + " " + R2 + \
                " --db " + nuc_ref_dir + \
                " --report " + base_id_out_fldr + base_id + "-nuc.kreport"

    # add flat flags
    command += " " + flat_flags

    # execute nuc k2 mapping
    process = subprocess.run(command.split(), stdout=subprocess.PIPE)
    if process.returncode != 0:
        raise Exception("Error running k2", command, "Error:", process.returncode)

    # Set up protein command 
    command = "kraken2" + \
                " --paired " + R1 + " " + R2 + \
                " --db " + prot_ref_dir + \
                " --report " + base_id_out_fldr + base_id + "-prot.kreport"

    # Need to correct for confidence score since K2 takes confidence as all 6 frames
    prot_flat_flats = ""
    adjust_confidence = False
    for flag in other_input_flags:
        if flag == "--confidence":
            adjust_confidence = True
        elif adjust_confidence == True:
            flag = str(float(flag) / 6.0)
            adjust_confidence = False
        prot_flat_flats += flag + " "
    if len(prot_flat_flats) > 0:
        prot_flat_flats = prot_flat_flats[:-1]

    # add protein flat flags
    command += " " + prot_flat_flats

    # execute prot k2 mapping
    process = subprocess.run(command.split(), stdout=subprocess.PIPE)
    if process.returncode != 0:
        raise Exception("Error running k2", command, "Error:", process.returncode)


# Queue up Bracken runs
for base_id in sampleid_fastq_dict:
    R1, R2 = sampleid_fastq_dict[base_id]

    # create sample-specific output folder
    base_id_out_fldr = out_dir + base_id + "/"
    os.makedirs(base_id_out_fldr, exist_ok=True)

    # Run bracken at both S and G levels
    command = bracken_exe + \
                " -d " + nuc_ref_dir + \
                " -i " + base_id_out_fldr + base_id + "-nuc.kreport" + \
                " -o " + base_id_out_fldr + base_id + "-S.bracken" + \
                " -r 150 -t 10" + \
                " -l " + "S"
    bracken_commands_to_run.append(p.apply_async(run_command, (command,)))

    command = bracken_exe + \
                " -d " + nuc_ref_dir + \
                " -i " + base_id_out_fldr + base_id + "-nuc.kreport" + \
                " -o " + base_id_out_fldr + base_id + "-G.bracken" + \
                " -r 150 -t 10" + \
                " -l " + "G"
    bracken_commands_to_run.append(p.apply_async(run_command, (command,)))

# Run Bracken in parallel process
# Initiate multiprocessing
[r.get() for r in bracken_commands_to_run]

# calculate alpha diversity from bracken results
for base_id in sampleid_fastq_dict:
    base_id_out_fldr = out_dir + base_id + "/"
    cal_alpha_diversity(base_id_out_fldr, base_id)

# Convert both K2 and Bracken outputs to feature tables across all samples
if restrict_taxids_to_mb == "TRUE":
    # Only spend time loading taxtree if we want to restrict TaxIDs in output
    print("\n\nLoading Taxtree")
    import taxtree
    tree = taxtree.TaxonomyTree(taxtree_loc)
    print("taxtree loaded\n\n")        
else:
    tree = ""

create_feature_table(sampleid_fastq_dict, out_dir, restrict_taxids_to_mb, tree)
