#!/usr/bin/env python3

import os

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
