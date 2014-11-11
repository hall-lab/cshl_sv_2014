#!/bin/bash -e

# ===============================================
# Cold Spring Harbor Laboratory Course
# Advanced Sequencing Technologies & Applications
# November 11-23, 2014
# ===============================================

# -----------------------------------------------
# 0. View BAM in IGV

# - Navigate to known deletion, and observe the
#     alignment patterns around it.

# -----------------------------------------------
# 1. Examine library insert size distribution

samtools view -f 0x0023 -F 0x051c NA12878.20.bam \
    | sed -n '1000001,2000000p;2000000q' \
    | awk '{ if ($7=="=") print $9 }' \
    > insert_size.txt

Rscript scripts/plot_histogram.R insert_size.txt
open insert_size.txt.pdf

# - Draw vertical lines on plot to demonstrate threshold
#     for discordant pairs
# - Show a poor quality insert distribution (bimodal)
#     for comparison. Draw vertical lines on this to
#     demonstrate reduction in power

# -----------------------------------------------
# 2. Generate insert size histogram for LUMPY

samtools view NA12878.20.bam \
    | sed -n '1000001,2000000p;2000000q' \
    | scripts/pairend_distro.py \
    -r 101 \
    -X 4 \
    -N 1000000 \
    -o NA12878.20.histo

# -----------------------------------------------
# 3. Extract split reads from BAM

samtools view -h NA12878.20.bam \
    | scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > NA12878.20.splitters.bam

# -----------------------------------------------
# 4. Extract discordant read pairs from BAM

samtools view -u -F 0x0002 NA12878.20.bam \
    | samtools view -u -F 0x0100 - \
    | samtools view -u -F 0x0004 - \
    | samtools view -u -F 0x0008 - \
    | samtools view -b -F 0x0400 - \
    > NA12878.20.discordants.bam

# -----------------------------------------------
# 5. Run LUMPY

lumpy \
    -mw 4 \
    -tt 0 \
    -pe bam_file:NA12878.20.discordants.bam,histo_file:NA12878.20.histo,mean:318.977050061,stdev:73.7166132046,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,id:10,min_mapping_threshold:20 \
    -sr bam_file:NA12878.20.splitters.bam,back_distance:10,min_mapping_threshold:20,weight:1,id:11,min_clip:20 \
    > NA12878.20.sv.bedpe

# -----------------------------------------------
# 6. Filtering SV calls

# - Evidence type (PE, SR, PE+SR)
# - Support
# - Overlap with 1000 Genomes calls
# - Long-read validation?
# - SVs that affect a gene

# -----------------------------------------------
# 7. Format conversions

# - Format for IGV
# - Convert to VCF

# -----------------------------------------------
# 8. Read-depth? Genotyping?

# -----------------------------------------------
# 9. Examine a pre-generated BEDPE of the entire
#     NA12878 genome?

# -----------------------------------------------
# 10. Examine a pre-generated BEDPE of a
#     tumor-normal pair?

# - Filter for somatic SVs affecting a gene








