#!/bin/bash

# Further optimization opportunities:
# Fastp maxes out at 16 threads -> multiple processes?
# Bracken is single threaded -> multiple processes?

# Define some variables/paths
threads=40
job_ID="run5_12.20.23"

trimmomatic_str="ILLUMINACLIP:/adapters/NexteraPE-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 SLIDINGWINDOW:5:20 MINLEN:50"
k2_hostr_conf=0.2
k2_class_conf=0.1

host_db="/dev/shm/k2_indexes/mouse_human/" 
nuc_ref_db="/dev/shm/k2_indexes/nt/" 
prot_ref_db="/dev/shm/k2_indexes/uni100_domain_k2_index/"

fastq_loc="../fastqs/"
raw_loc=$fastq_loc"raw/"
trim_loc=$fastq_loc"trim/"
hostr_loc=$fastq_loc"hostr/"

taxdmp_loc="../taxdmp_10.23.23/"

# Tie results directory to a job ID
result_dir="../results_job-"$job_ID"/"

if ! [[ -d $result_dir ]]; then
  echo -e "\nCreating output directory: ${result_dir}\n" >&2 
  mkdir -p "${result_dir}"
 fi

qual_rpt_loc=$result_dir"qual_rpts/"
if ! [[ -d $qual_rpt_loc ]]; then
  mkdir -p "${qual_rpt_loc}"
 fi

raw_qual=$qual_rpt_loc"raw_qc/"
trim_qual=$qual_rpt_loc"trim_qc/"
hostr_qual=$qual_rpt_loc"hostr_qc/"

classification_loc=$result_dir"taxclass/"

summarization_loc=$result_dir$job_ID"-summary.csv"

# Write out a job log file to a bash file
printf "$job_ID
ResultDir:$result_dir

HostDB:$host_db
NucDB:$nuc_ref_db
ProtDB:$prot_ref_db

TaxDMP:$taxdmp_loc

Trimmomatic Parameters:$trimmomatic_str
K2 Hostr Confidence:$k2_hostr_conf
K2 Classification Confidence:$k2_class_conf"\
> "${result_dir}/job_parameters.txt"

# Run qual reports on all raw FASTQs
echo "Running QC on raw reads"
Nano-QualReport/run $raw_loc $raw_qual $threads

# Run read trimming
echo "Running read trimming"
# Trimmomatic
Nano-ReadTrim/run trimmomatic TRUE $raw_loc $trim_loc \
                        -threads $threads \
                        $trimmomatic_str

#Fastp
#Nano-ReadTrim/run TRUE $raw_loc $trim_loc \
#        --cut_front --cut_front_window_size 1 --cut_front_mean_quality 3 \
#        --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 3 \
#        --cut_right --cut_right_window_size 5 --cut_right_mean_quality 20 \
#        -y -l 50 \
#        --thread $threads

# Run qual reports on all trimmed FASTQs
echo "Running QC on trimmed reads"
Nano-QualReport/run $trim_loc $trim_qual $threads

# Run host removal
echo "Running host removal"
Nano-HostRemoval/run kraken2 TRUE $trim_loc $hostr_loc $host_db \
        --threads $threads \
        --confidence $k2_hostr_conf --memory-mapping

# Run qual reports on all host removed FASTQs
echo "Running quality reports on host removed reads"
Nano-QualReport/run $hostr_loc $hostr_qual $threads

# Run K2 classification & bracken
echo "Running K2/bracken classification"
Nano-K2TaxClass/run TRUE $taxdmp_loc $hostr_loc $classification_loc $nuc_ref_db $prot_ref_db \
        --threads $threads \
        --confidence $k2_class_conf --memory-mapping

# Run post-run summary of read counts
echo "Running post-run summary for read counts"
Nano-SummarizeUtils/run $qual_rpt_loc $classification_loc $summarization_loc
