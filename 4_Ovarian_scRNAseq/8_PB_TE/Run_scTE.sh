#!/bin/bash

################################################################################
# scTE Run Script for Aging, VCD, and Foxl2 models
################################################################################

INDEX="/Users/minhookim/Programs/scTE/Genome_index/mm10/mm10.exclusive.idx"

# Function to run scTE
run_scTE () {
    local bam_file="$1"
    local sample_name="$2"

    ${SC_TE_CMD} -i "${bam_file}" \
                 -o "./scTE_out_${sample_name}" \
                 -x "${INDEX}" \
                 --hdf5 True \
                 -CB CR \
                 -UMI UR
}

################################################################################
# 1. Aging model
################################################################################

run_scTE "/Volumes/OIProject_I/1_Cellranger/Benayoun_lab/AC/YF_1/outs/possorted_genome_bam.bam" "YF_1"
run_scTE "/Volumes/OIProject_I/1_Cellranger/Benayoun_lab/AC/YF_2/outs/possorted_genome_bam.bam" "YF_2"
run_scTE "/Volumes/OIProject_I/1_Cellranger/Benayoun_lab/AC/OF_1/outs/possorted_genome_bam.bam" "OF_1"
run_scTE "/Volumes/OIProject_I/1_Cellranger/Benayoun_lab/AC/OF_2/outs/possorted_genome_bam.bam" "OF_2"

################################################################################
# 2. VCD model 
################################################################################

VCD_DIR="/Volumes/OIProject_I/1_Cellranger/Benayoun_lab/VCD"
for folder in "${VCD_DIR}"/*/; do
    [ -d "$folder" ] || continue
    bam_file="${folder}/outs/possorted_genome_bam.bam"
    sample_name=$(basename "$folder")
    run_scTE "${bam_file}" "${sample_name}"
done

################################################################################
# 3. Foxl2 model
################################################################################

FOXL2_DIR="/Volumes/OIProject_I/1_Cellranger/Benayoun_lab/Foxl2"
for folder in "${FOXL2_DIR}"/*/; do
    [ -d "$folder" ] || continue
    bam_file="${folder}/outs/possorted_genome_bam.bam"
    sample_name=$(basename "$folder")
    run_scTE "${bam_file}" "${sample_name}"
done
