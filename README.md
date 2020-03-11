# Project Description

A brief description of what this repository is for and what it contains

# Contributors
Michiel Smit smit2@bu.edu
Simran Makwana makwana@bu.edu
Emily Hughes ehug@bu.edu
Sumiti Sandhu sandhu08@bu.edu

# Repository Contents
run_extract.qsub: extract short reads to FASTQ format

quality_control.sh: process FASTQ files and extract quality metrics

run_tophat.qsub: utilize the Tophat tool to map splice junctions to the reads. This file mapped the reads from the P0 data to the mm9 mouse reference genome. Must have the following modules loaded: samtools, bowtie, and tophat

run_genebody.qsub: utilize the geneBody_coverage.py tool from RSeQC module. The read coverage over the gene body for the P0 data was calculated. Must have the following modules loaded: python3, samtools, and rseqc

run_cufflinks.qsub: utilize the Cufflinks tool to conduct differential expression analysis. This analysis was conducted on the P0 data. Must have the cufflinks module loaded

run_cuffdiff.qsub: utlize the Cuffdiff tool to conduct differential expression analysis. This was used to find differences in expression between P0 and Ad data. Must have the cufflinks module loaded
