#!/bin/bash

# Define some variables/paths
threads=8
job_ID="Fredrick_1"
result_dir="../results_job-"$job_ID"/"

trimmomatic_str="ILLUMINACLIP:/adapters/TruSeq2-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50"

# Read Locations
fastq_loc="../sequence_reads/Fredrick_1/"
raw_loc=$fastq_loc"raw/"
trim_loc=$fastq_loc"trim/"

# Create result directory if it doesn't already exist
if ! [[ -d $result_dir ]]; then
  echo -e "\nCreating output directory: ${result_dir}\n" >&2 
  mkdir -p "${result_dir}"
 fi

# Create quality report directories (FastQC, MultiQC)
qual_rpt_loc=$result_dir"qual_rpts/"
if ! [[ -d $qual_rpt_loc ]]; then
  mkdir -p "${qual_rpt_loc}"
 fi

raw_qual=$qual_rpt_loc"raw_qc/"
trim_qual=$qual_rpt_loc"trim_qc/"

# Define assembly output location
assembly_loc=$result_dir"assembly_out/"
if ! [[ -d $qual_rpt_loc ]]; then
  mkdir -p "${qual_rpt_loc}"
 fi

# Write out a job log file to a bash file
printf "$job_ID
ResultDir:$result_dir

Trimmomatic Parameters:$trimmomatic_str"
> "${result_dir}/job_parameters.txt"

"""
# Run qual reports on all raw FASTQs
echo "Running QC on raw reads"
Nano-QualReport/run $raw_loc $raw_qual $threads

# Run read trimming
echo "Running read trimming"
# Trimmomatic
Nano-ReadTrim/run trimmomatic TRUE $raw_loc $trim_loc \
                        -threads $threads \
                        $trimmomatic_str

# Run qual reports on all trimmed FASTQs
echo "Running QC on trimmed reads"
Nano-QualReport/run $trim_loc $trim_qual $threads
"""

# Conduct assembly
echo "Running assembly"
Nano-Assembly/run $trim_loc $assembly_loc -t $threads --rna
