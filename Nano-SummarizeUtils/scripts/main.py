#!/usr/bin/env python3

import sys
import subprocess
import os
import csv
from gen_plots import gen_count_boxplot, gen_domain_count_barplot

fastqc_parent_fldr = "/qual_rpts/"
taxclass_parent_fldr = "/taxclass/"
result_fldr = "/results/"

result_file = sys.argv[1]
result_loc = result_fldr + result_file
alpha_diversity_loc = result_fldr + result_file.split(".csv")[0] + "-alpha_diversity.csv"
boxplot_loc = result_fldr + result_file.split(".csv")[0] + "-counts.png"
domain_barplot_loc = result_fldr + result_file.split(".csv")[0] + "-domains.png"

master_readcount_dict = {}
desired_readcount_type_order = ["raw", "trim", "hostr", \
                                "k2_prot", "k2_nuc", "bracken_S", "bracken_G", \
                                "k2_prot_contains_mammals", "k2_nuc_contains_mammals"]

# First, gather all available FastQC reports for read counts
for fldr in os.listdir(fastqc_parent_fldr):
    if "raw" in fldr:
        readcount_type = "raw"
    elif "trim" in fldr:
        readcount_type = "trim"
    elif "hostr" in fldr:
        readcount_type = "hostr"
    else:
        raise Exception("Encountered unrecognized qual report fldr:", fldr)

    general_stats_loc = fastqc_parent_fldr + fldr + "/multiqc_data/multiqc_general_stats.txt"
    if os.path.isfile(general_stats_loc) == False:
        raise Exception("Unable to find read counts. Expected file loc:", general_stats_loc)

    with open(general_stats_loc, "r") as f:
        next(f)
        for row in f:
            row = row.replace("\n","").split("\t")
            sample_id = row[0]
            count = float(row[6])

            # Only grab forward read for sample tracking
            if "_R1" in sample_id:
                sample_id = sample_id.split("_R1")[0]
            else:
                continue

            if readcount_type == "hostr":
                sample_id = sample_id.split("-trim-hostr")[0]
            elif readcount_type == "trim":
                sample_id = sample_id.split("-trim")[0]

            if sample_id not in master_readcount_dict:
                master_readcount_dict[sample_id] = {}
            master_readcount_dict[sample_id][readcount_type] = count
    f.close()

# Next, gather kraken2 classification counts 
# For domain-level plotting, get counts per domain TaxID
sample_domain_taxid_counts = {}
for fldr in os.listdir(taxclass_parent_fldr):
    sample_id = fldr.split("-trim")[0]

    for filename in os.listdir(taxclass_parent_fldr + fldr):
        if "-prot.kreport" in filename:
            readcount_type = "k2_prot"
            taxcheck_type = "k2_prot_contains_mammals"
        elif "-nuc.kreport" in filename:
            readcount_type = "k2_nuc"
            taxcheck_type = "k2_nuc_contains_mammals"
        elif "-S.bracken" in filename:
            readcount_type = "bracken_S"
            taxcheck_type = ""
        elif "-G.bracken" in filename:
            readcount_type = "bracken_G"
            taxcheck_type = ""
        else:
            # Note - not else raise exception to deal w/ extra bracken report files
            continue

        # initialize sample_id in tracker
        if sample_id not in master_readcount_dict:
            master_readcount_dict[sample_id] = {}

        # kraken2 report format
        if "k2" in readcount_type:
            # initialize domain-level tracking if k2_prot
            if readcount_type == "k2_prot":
                sample_domain_taxid_counts[sample_id] = {}

            # extract root # of reads (this is # classified)
            num_classified = ""
            contains_mammals = False
            with open(taxclass_parent_fldr + fldr + "/" + filename, "r") as f:
                for row in f:
                    row = row.replace("\n","").split("\t")
                    perc, num_clade, num_taxid, rank, taxid, name = row

                    if name == "root":
                        num_classified = float(num_clade)
                    elif name == "Mammalia":
                        contains_mammals = True
                    
                    if rank == "D" and readcount_type == "k2_prot":
                        sample_domain_taxid_counts[sample_id][name] = float(num_clade)
            f.close()

            if num_classified == "":
                raise Exception("Unable to extract # of classified reads for report:", taxclass_parent_fldr + fldr + "/" + filename)
        
            # Add mammal tracking for k2 rpts
            master_readcount_dict[sample_id][taxcheck_type] = contains_mammals

        # bracken report format
        if "bracken" in readcount_type:
            # sum total # of reads
            num_classified = 0
            with open(taxclass_parent_fldr + fldr + "/" + filename, "r") as f:
                next(f)
                for row in f:
                    row = row.replace("\n","").split("\t")
                    name, taxid, rank, k2_assigned, added, new_est, fract_total = row

                    num_classified += float(new_est)
            f.close()

        # Add value to the tracking dict
        master_readcount_dict[sample_id][readcount_type] = num_classified

# now, write the output out!
with open(result_loc, "w") as f:
    out = csv.writer(f)

    header_row = ["Sample ID"] + desired_readcount_type_order
    out.writerow(header_row)

    for sample_id in master_readcount_dict:
        sample_counts = master_readcount_dict[sample_id]

        out_row = [sample_id]
        for readcount_type in desired_readcount_type_order:
            if readcount_type in sample_counts:
                count = sample_counts[readcount_type]
            else:
                count = ""
            out_row.append(count)

        out.writerow(out_row)
f.close()

# Get diversity metrics 
alpha_types = ["BP", "Sh", "F", "Si", "ISi"]
sample_diversity_dict = {}

# Next, gather diversity metrics
for fldr in os.listdir(taxclass_parent_fldr):
    sample_id = fldr.split("-trim")[0]

    if sample_id == "summary":
        continue

    sample_diversity_dict[sample_id] = {}

    for filename in os.listdir(taxclass_parent_fldr + fldr):
        if "-alpha_diversity.csv" in filename:
            with open(taxclass_parent_fldr + fldr + "/" + filename, "r") as f:
                reader = csv.reader(f)
                next(reader, None)

                for row in reader:
                    alpha, diversity = row

                    if diversity != "":
                        diversity = float(diversity.replace('\\n"',''))

                        sample_diversity_dict[sample_id][alpha] = diversity
            f.close()
    
# write output file if diversity was calculated
if len(sample_diversity_dict) > 0:
    with open(alpha_diversity_loc, "w") as f:
        out = csv.writer(f)

        header_row = ["Sample_ID"]
        header_row = header_row + alpha_types
        out.writerow(header_row)

        for sample_id in sample_diversity_dict:
            sample_alphas = sample_diversity_dict[sample_id]

            out_row = [sample_id]
            for alpha in alpha_types:
                if alpha in sample_alphas:
                    value = sample_alphas[alpha]
                else:
                    value = ""
                out_row.append(value)
            
            out.writerow(out_row)
    f.close()

# Generate output count boxplot export
gen_count_boxplot(result_loc, boxplot_loc)

# Generate domain-level distributions barplot export
gen_domain_count_barplot(sample_domain_taxid_counts, domain_barplot_loc)
