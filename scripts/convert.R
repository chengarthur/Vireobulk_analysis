#!/usr/bin/env Rscript
##optparses
ConvertoPdFormat<-function(csvdf){
  df<-data.frame("variants"=paste(csvdf$Chr,csvdf$Start,csvdf$Ref,csvdf$Alt,sep = "_"),"function"=csvdf$Func.refGene,"genes"=csvdf$Gene.refGene,"avsnp"=csvdf$avsnp150)
  return(df)
}


library(optparse,quietly = T)
parser = OptionParser(description='R script for carculate biotype')
parser = add_option(parser, '--annofile', default = 'filtered_table.hg38_multianno.csv', type = "character", help = "vcf file annotated by annovar")
parser = add_option(parser, '--out',default = 'scvcf_annotated_processed.csv', type = "character", help = "out file name ")
args <- parse_args(parser)
df=args$annofile
testdf<-read.csv(file =df)
write.csv(ConvertoPdFormat(testdf),args$out)


