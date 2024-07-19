#!/usr/bin/env python3

import os
import csv
import subprocess


def cal_alpha_diversity(base_id_out_fldr, base_id):
    """
    """
    alpha_div_exe = "/KrakenTools/DiversityTools/alpha_diversity.py"

    # For the bracken output, calculate alpha diversity measurements
    diversity_out = base_id_out_fldr + base_id + "-alpha_diversity.csv"

    alpha_types = ["BP", "Sh", "F", "Si", "ISi"]

    # get species-level a diversities
    bracken_out_file = base_id_out_fldr + base_id + "-S.bracken"

    with open(diversity_out, "w") as f:
        out = csv.writer(f)
        out.writerow(["Alpha Diversity Metric", "Value"])

        for alpha in alpha_types:
            command = alpha_div_exe + \
                        " -f " + bracken_out_file + \
                        " -a " + alpha
            process = subprocess.run(command.split(), stdout=subprocess.PIPE)

            subprocess_out = str(process.stdout).split("\n")
            # the diversity follows the :
            diversity = ""
            for row in subprocess_out:
                if ": " in row:
                    diversity = row.split(": ")[1]
            
            if diversity == "":
                print("Unable to compute diversity for:", base_id, " diversity type:", alpha)

            out.writerow([alpha, diversity])
    f.close()

    return

def run_command(command):
    #process = subprocess.run(command.split(), stdout=subprocess.PIPE)
    #if process.returncode != 0:
    #    print("Error running bracken", command, "Error:", process.returncode)
    os.system(command)
    return


def search_for_threads_in_flags(other_input_flags):
    """
    """
    find_threads = False
    num_threads = 1
    for flag in other_input_flags:
        if flag == "--threads":
            find_threads = True
        elif find_threads == True:
            num_threads = int(flag)
            find_threads = False

    return num_threads



def find_sample_fastqs(fastq_dir):
    """
    """
    sampleid_fastq_dict = {}

    for filename in os.listdir(fastq_dir):
        if "_R1" in filename:
            base_id = filename.split("_R1")[0]
            index = 0
        elif "_R2" in filename:
            base_id = filename.split("_R2")[0]
            index = 1
        else:
            print("WARNING - Encountered unacceptable filetype for PE reads:", filename)

        if base_id not in sampleid_fastq_dict:
            sampleid_fastq_dict[base_id] = ["", ""]

        sampleid_fastq_dict[base_id][index] = fastq_dir + filename

    # check sampleid_fastq_dict for base_ids without 2 FASTQs
    for base_id in sampleid_fastq_dict:
        R1, R2 = sampleid_fastq_dict[base_id]

        if R1 == "" or R2 == "":
            raise Exception("Did not find both R1 and R2 for base_id:", base_id)

    return sampleid_fastq_dict


def check_if_taxid_in_desired_clade(og_taxid, tree):
    """
    """
    desired_tax_clades = [
        2, \
        4751, \
        10239
    ]

    if int(og_taxid) not in tree:
        print("TaxID not in taxid:", og_taxid)
        return False

    lineage = tree.ascend(int(og_taxid))
    lineage_has_desired_clade = False

    for node in lineage:
        taxid = node.taxid

        if taxid in desired_tax_clades:
            lineage_has_desired_clade = True

    return lineage_has_desired_clade


def load_taxid_lists_bracken(sample_kreports, restrict_taxids_to_mb, tree, \
                                out_dir, bracken_type):
    """
    """
    taxid_name_dict = {}
    base_id_counts = {}
    all_uniq_taxids = set()
    sample_taxid_abund_tracker = {}

    for base_id, kreport_loc in sample_kreports:
        sample_taxid_tracker = []
        with open(kreport_loc, "r") as f:
            next(f)
            for row in f:
                row = row.replace("\n","").split("\t")

                name, taxid, taxlvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads = row
                sample_taxid_tracker.append([name, taxid, taxlvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads])
        f.close()

        count_classified = 0
        desired_taxid_counts = []
        for name, taxid, taxlvl, kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads in sample_taxid_tracker:
            if restrict_taxids_to_mb == "TRUE":
                lineage_has_desired_clade = check_if_taxid_in_desired_clade(taxid, tree)
            else:
                lineage_has_desired_clade = True
            
            if lineage_has_desired_clade == False:
                continue
            
            count_classified += int(new_est_reads)
            desired_taxid_counts.append([taxid, int(new_est_reads)])
            taxid_name_dict[taxid] = name
        
        base_id_counts[base_id] = count_classified
        sample_taxid_abund_tracker[base_id] = {}
        for taxid, count in desired_taxid_counts:
            fract = count / count_classified
            sample_taxid_abund_tracker[base_id][taxid] = fract
            all_uniq_taxids.add(taxid)

    # Now go through samples, creating 2D sample-RA matrix 
    # Fill in RA of 0 if taxid not observed in sample
    all_uniq_taxids = list(all_uniq_taxids)
    final_feature_dict = {}
    for sample_id in sample_taxid_abund_tracker:
        final_feature_dict[sample_id] = []

        sampleid_taxid_dict = sample_taxid_abund_tracker[sample_id]

        for taxid in all_uniq_taxids:
            if taxid in sampleid_taxid_dict:
                fract = sampleid_taxid_dict[taxid]
            else:
                fract = 0.0

            final_feature_dict[sample_id].append(fract)

    # Write out feature table
    with open(out_dir + bracken_type + "_" "bracken-feature-table.csv", "w") as f:
        out = csv.writer(f)

        # header
        header_row = ["Sample_ID"]
        for taxid in all_uniq_taxids:
            header_row.append(taxid)
        out.writerow(header_row)

        for sample_id in sample_taxid_abund_tracker:
            taxid_array = final_feature_dict[sample_id]

            out_row = [sample_id] + taxid_array
            out.writerow(out_row)
    f.close()

    return base_id_counts


def load_taxid_lists_kraken(sample_kreports, restrict_taxids_to_mb, tree, \
                                out_dir, k2_type):
    """
    """
    taxid_name_dict = {}
    base_id_counts = {}
    all_uniq_taxids = set()
    sample_taxid_abund_tracker = {}

    for base_id, kreport_loc in sample_kreports:
        sample_taxid_tracker = []
        with open(kreport_loc, "r") as f:
            for row in f:
                row = row.replace("\n","").split("\t")

                perc, num_clade, num_taxid, rank, taxid, name = row
                sample_taxid_tracker.append([perc, num_clade, num_taxid, rank, taxid, name])
        f.close()

        total_sample_count = 0
        count_classified = 0
        desired_taxid_counts = []
        for perc, num_clade, num_taxid, rank, taxid, name in sample_taxid_tracker:
            total_sample_count += int(num_taxid)

            if rank != "S":
                continue

            if restrict_taxids_to_mb == "TRUE":
                lineage_has_desired_clade = check_if_taxid_in_desired_clade(taxid, tree)
            else:
                lineage_has_desired_clade = True
            
            if lineage_has_desired_clade == False:
                continue
            
            count_classified += int(num_clade)
            desired_taxid_counts.append([taxid, int(num_clade)])
            taxid_name_dict[taxid] = name
        
        base_id_counts[base_id] = [count_classified, total_sample_count]
        sample_taxid_abund_tracker[base_id] = {}
        for taxid, count in desired_taxid_counts:
            fract = count / count_classified
            sample_taxid_abund_tracker[base_id][taxid] = fract
            all_uniq_taxids.add(taxid)

    # Now go through samples, creating 2D sample-RA matrix 
    # Fill in RA of 0 if taxid not observed in sample
    all_uniq_taxids = list(all_uniq_taxids)
    final_feature_dict = {}
    for sample_id in sample_taxid_abund_tracker:
        final_feature_dict[sample_id] = []

        sampleid_taxid_dict = sample_taxid_abund_tracker[sample_id]

        for taxid in all_uniq_taxids:
            if taxid in sampleid_taxid_dict:
                fract = sampleid_taxid_dict[taxid]
            else:
                fract = 0.0

            final_feature_dict[sample_id].append(fract)

    # Write out feature table
    with open(out_dir + k2_type + "_" "k2-feature-table.csv", "w") as f:
        out = csv.writer(f)

        # header
        header_row = ["Sample_ID"]
        for taxid in all_uniq_taxids:
            header_row.append(taxid)
        out.writerow(header_row)

        for sample_id in sample_taxid_abund_tracker:
            taxid_array = final_feature_dict[sample_id]

            out_row = [sample_id] + taxid_array
            out.writerow(out_row)
    f.close()

    return base_id_counts


def create_feature_table(sampleid_fastq_dict, out_dir, restrict_taxids_to_mb, tree):
    """
    """
    base_id_nuc_kreports = []
    base_id_prot_kreports = []
    base_id_bracken_S_rpts = []
    base_id_bracken_G_rpts = []
    for base_id in sampleid_fastq_dict:
        base_id_out_fldr = out_dir + base_id + "/"

        nuc_k2_report = base_id_out_fldr + base_id + "-nuc.kreport"
        base_id_nuc_kreports.append([base_id, nuc_k2_report])

        prot_k2_report = base_id_out_fldr + base_id + "-prot.kreport"
        base_id_prot_kreports.append([base_id, prot_k2_report])

        # Note - if no classifications in kraken2 (e.g., blanks), bracken 
        # won't produce an output file. In this case, skip!
        bracken_S_rpt = base_id_out_fldr + base_id + "-S.bracken"
        if os.path.isfile(bracken_S_rpt) == True:
            base_id_bracken_S_rpts.append([base_id, bracken_S_rpt])

        bracken_G_rpt = base_id_out_fldr + base_id + "-G.bracken"
        if os.path.isfile(bracken_G_rpt) == True:
            base_id_bracken_G_rpts.append([base_id, bracken_G_rpt])

    # Assign a summary-level location
    summary_report_loc = out_dir + "summary/"
    os.makedirs(summary_report_loc, exist_ok=True)

    nuc_base_id_counts = load_taxid_lists_kraken(base_id_nuc_kreports, \
                            restrict_taxids_to_mb, tree, \
                            summary_report_loc, "nuc")
    prot_base_id_counts = load_taxid_lists_kraken(base_id_prot_kreports, \
                            restrict_taxids_to_mb, tree, \
                            summary_report_loc, "prot")
    if len(base_id_bracken_S_rpts) > 0:
        bracken_S_counts = load_taxid_lists_bracken(base_id_bracken_S_rpts, \
                                restrict_taxids_to_mb, tree, \
                                summary_report_loc, "S")
    else:
        bracken_S_counts = {}
        print("WARNING - no Bracken S reports")
    if len(base_id_bracken_G_rpts) > 0:
        bracken_G_counts = load_taxid_lists_bracken(base_id_bracken_G_rpts, \
                                restrict_taxids_to_mb, tree, \
                                summary_report_loc, "G")
    else:
        bracken_G_counts = {}
        print("WARNING - no Bracken G reports")
    
    with open(summary_report_loc + "classification-count-summary.csv", "w") as f:
        out = csv.writer(f)
        out.writerow(["Base ID", "Tot Final Seq Count", \
                            "Nuc K2 Classified Count", \
                            "Prot K2 Classified Count", \
                            "Bracken S Classified Count", \
                            "Bracken G Classified Count"])

        for base_id in nuc_base_id_counts:
            nuc_k2_count_classified, total_sample_count = nuc_base_id_counts[base_id]
            prot_k2_count_classified, _ = prot_base_id_counts[base_id]

            if base_id in bracken_S_counts:
                bracken_s_count = bracken_S_counts[base_id]
            else:
                bracken_s_count = 0
            if base_id in bracken_G_counts:
                bracken_g_count = bracken_G_counts[base_id]
            else:
                bracken_g_count = 0

            out.writerow([base_id, total_sample_count, \
                            nuc_k2_count_classified, \
                            prot_k2_count_classified, \
                            bracken_s_count, bracken_g_count])
    f.close()

    return
