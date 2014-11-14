cd ~/HALL_2014/results

#----------------------------------------------------------------------------------------------------------
# 1) GENERATE THE INSERT SIZE DISTRIBUTION FILE
#----------------------------------------------------------------------------------------------------------

# This is what LUMPY uses to construct the lumpy probability distributions used for clustering
# -r is the read length, -N is the number of alignments to use
# -X is the the number of standard deviations, -o is the output file
samtools view NA12878.20.bam \
    | head -n 1000000 \
    | tail -n 100000 \
    | ~/HALL_2014/bin/pairend_distro.py -r 101 -X 4 -N 100000 -o histo.out

# If you want to be extra careful, plot LUMPY's determination of the size distribution and compare to our own.
# Notice this is smaller because it does not include one read length

# open a new terminal window (command + N) and open the R program by typing "r" in the terminal
setwd("~/HALL_2014/results")
histofile = read.table("histo.out")
dev.new(); plot(histofile[,1],histofile[,2], type="h")

#----------------------------------------------------------------------------------------------------------
# 2) EXTRACT DISCORDANT READ-PAIR MAPPINGS
#----------------------------------------------------------------------------------------------------------
# QUESTION: what does each line of this command do?
# ANSWER: not concordant; query not unmapped; mate not unmapped; not secondary; not a duplicate.
# QUESTION: what does -u and -b do?
# ANSWER: -u says to write uncompressed BAM, which make things faster because each line does not need to be
# uncompressed and compressed again before the next command in the stream.
# -b says to write compressed BAM; the default is to write SAM
samtools view -u -F 0x0002 NA12878.20.bam \
    | samtools view -u -F 0x4 - \
    | samtools view -u -F 0x8 - \
    | samtools view -u -F 0x100 - \
    | samtools view -b -F 0x400 - \
    > NA12878.20.discordants.bam &

#----------------------------------------------------------------------------------------------------------
# 3) EXTRACT SPLIT-READ MAPPINGS
#----------------------------------------------------------------------------------------------------------

# Start this command immediately since it takes ~5 minutes to run
samtools view -h NA12878.20.bam \
    | ~/HALL_2014/bin/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > NA12878.20.splitters.bam &

# Questions to talk about while the above command is running:
# a) How are split-read alignments reported in the bam file?
# b) What does extractSplitReads_BwaMem do?
# c) What is a soft-clipped vs. hard-clipped read?
# d) What fraction of aligned reads have split-read alignments?
    # There are 31330851 aligned reads in the original bam file that are primary and not duplicates
    # There are 158519 split-read alignments based on samtools flagstat =  0.0051 = 0.51%
# e) What fraction of readpairs are discordant alignments?
    # There are 15622923 non-duplicate readpairs with both reads mapped, and 341974 discordant readpairs 
    # = 0.0219 = 2.19%
# These values are typical. Datasets with a high frequency of discordant or split-read mappings have a problem.

#----------------------------------------------------------------------------------------------------------
# 4) RUN LUMPY AT TWO DIFFERENT PARAMETER SETTINGS
#----------------------------------------------------------------------------------------------------------

# First, run Lumpy with "naive" parameters, which make the assumption that SV detection is easy. 
# What does each parameter do? 
# Type  "~/HALL_2014/bin/lumpy -h" and hit return to see list of options
~/HALL_2014/bin/lumpy \
    -mw 2 \
    -tt 0 \
    -pe bam_file:NA12878.20.discordants.bam,histo_file:histo.out,mean:319.551326228,stdev:74.2952533362,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,id:10,min_mapping_threshold:20 \
    -sr bam_file:NA12878.20.splitters.bam,back_distance:10,min_mapping_threshold:20,weight:1,id:11,min_clip:20 \
    > naive.out

# Second, run Lumpy with "strict" parameters to map high-confidence breakpoints
# QUESTION: what are the two key differences between the naive and strict parameters that we are using? 
~/HALL_2014/bin/lumpy \
    -mw 7 \
    -tt 0 \
    -x ../annotations/exclude.b37.bed \
    -pe bam_file:NA12878.20.discordants.bam,histo_file:histo.out,mean:319.551326228,stdev:74.2952533362,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,id:10,min_mapping_threshold:20 \
    -sr bam_file:NA12878.20.splitters.bam,back_distance:10,min_mapping_threshold:20,weight:1,id:11,min_clip:20 \
    > strict.out
    
# How many SV breakpoint calls are in the two datasets?
wc -l naive.out strict.out
# 3380 naive.out
#  126 strict.out
# 3506 total

#?
# FYI, running the strict version without the exclude list yields 254
# FYI, running the naive version with the exclude list yields 2828 calls

# Reformat the raw lumpy output file to something more user friendly 
# (this will not be necessary with next LUMPY version, if Ryan Layer doesn't continue to ignore me)
# This requires a config file describing the sample IDs 
~/HALL_2014/bin/lumpyToBedpe -b naive.out -c ~/HALL_2014/bin/config.txt > breakpoints.naive.bedpe
~/HALL_2014/bin/lumpyToBedpe -b strict.out -c ~/HALL_2014/bin/config.txt > breakpoints.strict.bedpe

# Take a look at these files using "less"

# LUMPY bedpe file format:
# 1 = chromosome a; 
# 2 = coordinate start a (leftmost position of the first breakpoint-containing genomic interval) 
# 3 = coordinate end a (rightmost position of the first breakpoint-containing genomic interval) 
# 4 = chromosome b; 
# 5 = coordinate start b (leftmost position of the second breakpoint-containing genomic interval) 
# 6 = coordinate end b (rightmost position of the second breakpoint-containing genomic interval) 
# 7 = breakpoint ID; 
# 8 = support (total number of read-pair + split-read measurements)
# 9 = strand a (direction of breakpoint relative to read mappings; "+" indicates to right, "-" to left) 
# 10 = strand b (direction of breakpoint relative to read mappings; "+" indicates to right, "-" to left)
# 11 = variant type (DEL=deletion; DUP=duplication; INV=inversion; INT=inter-chromosomal) 
# 12 = evidence types detected (PE/SR)
# 13 = strand combinations clustered (useful for inversions and reciprocal translocations)
# 14 = sampleList (list of samples that have the breakpoint) 
# 15 = evidenceSampleList (detailed list of samples, evidence types, and evidence observations)
# 16-N: sample_N total support


# If you want to look at breakpoint calls for the full genome,
# use this file: ~/HALL_2014/supplemental/NA12878.lumpy.bedpe
# Note, however, that they were generated with the SpeedSeq pipeline, 
# which is slightly different, so there are slightly fewer calls (filtered false positives)

# QUESTION: how many of each variant type were detected in the naive and strict callsets?
cat breakpoints.naive.bedpe | cut -f 11 | sort | uniq -c
#  467 DEL
#   86 DUP
# 2390 INT
#  437 INV

cat breakpoints.strict.bedpe | cut -f 11 | sort | uniq -c
#  82 DEL
#  12 DUP
#  27 INT
#   5 INV

# Why are the two datasets so different?
