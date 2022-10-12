#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
from base import *
from vireoSNP.utils.vcf_utils import load_VCF, write_VCF, parse_donor_GPb
from version import __version__
import pandas as pd
import sys

def main():
   
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
   bulk_data=args.cell_data
   donor_data=args.donor_file
   mode_tag=args.mode
   
   
   print(bulk_data)
   print(donor_data)
   print(mode_tag)
   print(__version__)
   
   ##test1 
   
def preprocess(bulk_data,donor_data):
   bulk_data=load_VCF("/home/flyflyzhizhi/Vireobulk_analysis/data/cellSNP.base.vcf.gz", biallelic_only=True,load_sample=False,format_list=None,)
   donor_vcf=load_VCF(   "/home/flyflyzhizhi/Vireobulk_analysis/data/filter_pbmc10donors.vcf.gz",  biallelic_only=True,    load_sample=True, sparse=False , format_list=None)
   
 #  if model== "bulk":
  #     df_anno=pd.read_csv("data/scvcf_annotated_processed.csv")
   Ndonors=10
   df_b,df_d,m_gt=rematch_merge_vcf_noannotated(bulk_data,donor_vcf)
   df_b,df_d,m_gt=removedfault(df_b, df_d,m_gt)
   a=np.array(df_b["AD"])
   d=np.array(df_b["DP"])
   a=np.array(list(map(float,a)))
   d=np.array(list(map(float,d)))
   model =vireoSNP.VireoBulk(10)

   model.fit(a,d,m_gt,learn_theta=True )
   return (model)
   
def preprocess_with_anno():
    bulk_data=load_VCF("/home/flyflyzhizhi/Vireobulk_analysis/data/cellSNP.base.vcf.gz", biallelic_only=True,load_sample=False,format_list=None,)
    donor_vcf=load_VCF(   "/home/flyflyzhizhi/Vireobulk_analysis/data/filter_pbmc10donors.vcf.gz",  biallelic_only=True,   load_sample=True, sparse=False , format_list=None)
    df_anno=pd.read_csv("data/scvcf_annotated_processed.csv")
    Ndonors=10
    df_b,df_d,m_GT=rematch_merge_vcf_annotated(bulk_data,donor_vcf,df_anno)
    df_b,df_d,m_GT=removedfault(df_b,df_d,m_GT)
    a=np.array(df_b["AD"])
    d=np.array(df_b["DP"])
    a=np.array(list(map(float,a)))
    d=np.array(list(map(float,d)))
    model =vireoSNP.VireoBulk(10)
    gene_unique=np.unique(df_d["genes"].tolist())
    model.fit(a,d,m_GT,learn_theta=True )
    premodel=model.psi
    bulk_demuti_result= prediction_bulk(df_b,df_d,m_GT,donor_vcf,premodel,gene_unique)
    return(0)
    
main()
   
   
    
