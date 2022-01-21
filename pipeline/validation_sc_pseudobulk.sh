#!/bin/sh
#require：
#https://github.com/10XGenomics/subset-bam
#https://github.com/single-cell-genetics/cellsnp-lite
#whitelist_raw_sc_bam
subset-bam_linux --bam "cellranger_genome_bam.bam" --cell-barcodes whitelistbarcode.tsv  --log-level info --cores 20 --out-bam only_good_cell.bam

#random sampling & index
samtools index only_good_cell.bam
samtools view only_good_cell.bam -s 0.8 -b -o 0.8coverage_good_cell.bam 

#pileup
cellsnp-lite -s "only_good_cell.bam" -O pbmc10donor -R ref_donors.vcf.gz -b whitelistbarcode.tsv --gzip -p 15 &

