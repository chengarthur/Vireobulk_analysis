#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
#from .utils.base import VCF_to_sdf
from base import *
from vireoSNP.utils.vcf_utils import load_VCF, write_VCF, parse_donor_GPb
from version import __version__
import pandas as pd
import sys

def main():
   """
    parser = argparse.ArgumentParser(description='input bulk pile up result reference,and annotation for SNPs(optional)')
   
   
   ##sample options  data
   parser.add_argument("--cellData", "-c", dest="cell_data", default=None,required=True,
       help=("The cell genotype file in VCF format or cellSNP folder with "
             "sparse matrices."))
   
   parser.add_argument("--donorFile", "-d", dest="donor_file", default=None,required=True,
       help=("The donor genotype file in VCF format. Please filter the sample "
       "and region with bcftools -s and -R first!"))
   
   
   parser.add_argument("--outDir", "-o", dest="out_dir", default=None,
       help=("Dirtectory for output files [default: $cellFilePath/vireo]"))

   parser.add_argument("--annotation", "-a", dest="anno_file", default=None,
    help=("annotation for SNPs and Genes"))
    
   parser.add_argument("--forceLearnGT", dest="force_learnGT", default=False,
       action="store_true", help="If use, treat donor GT as prior only.")

   parser.add_argument('--gene', dest='mode', action='store_const',
                    const="gene", default="bulk",help='gene level demutiplexing (default: demutiplex at propotion level)')

   args = parser.parse_args()
   """
   print(__version__)
   ##test1 
   
#def preprocess():
   donor_vcf=load_VCF(   "/home/flyflyzhizhi/Vireobulk_analysis/data/filter_pbmc10donors.vcf.gz",  biallelic_only=True,    load_sample=True, sparse=False , format_list=None)
   bulk_data=load_VCF("/home/flyflyzhizhi/Vireobulk_analysis/data/cellSNP.base.vcf.gz", biallelic_only=True,load_sample=False,format_list=None,)
 #  if model== "bulk":
  #     df_anno=pd.read_csv("data/scvcf_annotated_processed.csv")
   df_b=merged_VCF_to_sdf(bulk_data)
   Ndonors=10
   df_d=merged_VCF_to_sdf(donor_vcf)
   donor_tensor = vireoSNP.vcf.parse_donor_GPb(donor_vcf['GenoINFO']["GT"], "GT")
   df_b,df_d,m_gt=removedfault(df_b, df_d, donor_tensor)
main()
   
   
    
