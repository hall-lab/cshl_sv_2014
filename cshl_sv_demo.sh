#!/bin/bash -e

# ===============================================
# Cold Spring Harbor Laboratory Course
# Advanced Sequencing Technologies & Applications
# November 11-23, 2014
# ===============================================

# -----------------------------------------------
# 1. Examine library insert size distribution

samtools view -f 0x0023 -F 0x051c NA12878.20.bam \
    | sed -n '1000000,2000000p;2000000q' \
    | awk '{ if ($7=="=") print $9 }' \
    > insert_size.txt

#
Rscript scripts/plot_histogram.R insert_size.txt
open insert_size.txt.pdf

# -----------------------------------------------
# 2. Generate insert size histograms for LUMPY

samtools view NA12878.20.bam \
    | sed -n '1000000,2000000p;2000000q' \
    | python scripts/pairend_distro.py \
    -r 101 \
    -X 4 \
    -N 1000000 \
    -o NA12878.20.histo

# -----------------------------------------------
# 2. Extract split reads from BAM

samtools view -h NA12878.20.bam \
    | scripts/extractSplitReads_BwaMem -i stdin \
    | awk '$0~"^@" || $7=="="' \
    | samtools view -Sb - \
    > NA12878.20.splitters.bam

# -----------------------------------------------
# 3. Extract discordant read pairs from BAM

samtools view -u -F 0x0002 NA12878.20.bam \
    | samtools view -u -F 0x0100 - \
    | samtools view -u -F 0x0004 - \
    | samtools view -u -F 0x0008 - \
    | samtools view -h -F 0x0400 - \
    | awk '$0~"^@" || $7=="="' \
    | samtools view -Sb - \
    > NA12878.20.discordants.bam

# -----------------------------------------------
# 5. Run LUMPY

lumpy \
    -mw 4 \
    -tt 0 \
    -pe bam_file:NA12878.20.discordants.bam,histo_file:NA12878.20.histo,mean:317.69018183,stdev:72.8912619975,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,id:10,min_mapping_threshold:20 \
    -sr bam_file:NA12878.20.splitters.bam,back_distance:10,min_mapping_threshold:20,weight:1,id:11,min_clip:20 \
    > NA12878.20.sv.bedpe
