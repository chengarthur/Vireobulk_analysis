#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse
from base import *
from vireoSNP.utils.vcf_utils import load_VCF, write_VCF, parse_donor_GPb
from version import __version__
import pandas as pd
import os

def main():
   
   parser = argparse.ArgumentParser(description='input bulk pile up result reference,and annotation for SNPs(optional)')
   
   
   ##sample options  data
   parser.add_argument("--cellData", "-c", dest="cell_data", default=None,required=True,
       help=("The cell genotype file in VCF format or cellSNP folder with "
             "sparse matrices."))
   
   parser.add_argument("--donorFile", "-d", dest="donor_file", default=None,required=True,
       help=("The donor genotype file in VCF format. Please filter the sample "
       "and region with bcftools -s and -R first! and the number of samples should be matched!"))
   
   
   parser.add_argument("--outDir", "-o", dest="out_dir", default=os.getcwd(),
       help=("Dirtectory for output files [default use os.getcwd as directory"))
   
   parser.add_argument("--preflixname", "-n", dest="pre_name", default=None,
       help=("preflix name of output names"))
   

   parser.add_argument("--annotation", "-a", dest="anno_file", default=None,
    help=("annotation for SNPs and Genes"))
    
   parser.add_argument("--notLearnGT", dest="force_learnGT", default=True,
       action="store_false", help="If use, close learnGT mode for ASE in donor ratio demutiplex, default use learned mode")
   ### when this flag is input , run gene level demutiplexing mode
   parser.add_argument('--gene', dest='mode', action='store_const',
                    const="gene", default="bulk",help='If use,into gene level demutiplexing (default: demutiplex at propotion level)')
##init args
   args = parser.parse_args()
   prename=args.pre_name
   bulk_data=args.cell_data
   donor_data=args.donor_file
   mode_tag=args.mode
   forceLearnGT_tag=args.force_learnGT
   outdir=args.out_dir
   anno_data=args.anno_file
   df=load_VCF( donor_data,biallelic_only=True,    load_sample=True, sparse=False , format_list=None)
   
   """breakpoint test
   print(bulk_data)
   print(donor_data)
   print(mode_tag)
   print(forceLearnGT_tag)
   print(outdir)
   print(__version__)
   """
   
   ##-c "/home/flyflyzhizhi/Vireobulk_analysis/data/cellSNP.base.vcf.gz" -d "/home/flyflyzhizhi/Vireobulk_analysis/data/filter_pbmc10donors.vcf.gz" 
   ## run bulk demutiplexing first
   if (mode_tag=="bulk"):
       model_returned=preprocess(bulk_data, donor_data,forceLearnGT_tag)
       print(model_returned.psi)
       OUTPUT_MODE=1
       modename=os.path.join(outdir,"%smodelsummary.txt") %prename
       with open(modename,"w") as f:
           f.write(','.join(df['samples']))
           f.write("\n")
           f.write(str(model_returned.psi).replace("\n",""))
           f.write("\n")
           f.write(str(model_returned.theta))
       
      ##-c "/home/flyflyzhizhi/Vireobulk_analysis/data/cellSNP.base.vcf.gz" -d "/home/flyflyzhizhi/Vireobulk_analysis/data/filter_pbmc10donors.vcf.gz" --gene -a "/home/flyflyzhizhi/Vireobulk_analysis/data/scvcf_annotated_processed.csv" 

   if (mode_tag=="gene"):
       model_returned,bulk_demuti_result=preprocess_with_anno(bulk_data, donor_data, anno_data,forceLearnGT_tag)
       print(model_returned.psi)
       print(bulk_demuti_result)
       OUTPUT_MODE=2
       modename=os.path.join(outdir,"%smodelsummary.txt")%prename
       with open(modename,"w") as f:
           f.write(','.join(df['samples']))
           f.write("\n")
           f.write(str(model_returned.psi).replace("\n",""))
           f.write("\n")
           f.write(str(model_returned.theta))
    
       bulk_demuti_result.to_csv(os.path.join(outdir,"%sdemutiplxing_result_output.csv"))%prename
       
   ##test1 
   
def preprocess(bulk_data1,donor_data1,forceLearnGT_tag):
    bulk_data=load_VCF(bulk_data1, biallelic_only=True,load_sample=False,format_list=None,)
    donor_vcf=load_VCF( donor_data1,  biallelic_only=True,    load_sample=True, sparse=False , format_list=None)
   

    Ndonors=len(donor_vcf['samples'])
    df_b,df_d,m_gt=rematch_merge_vcf_noannotated(bulk_data,donor_vcf)
    df_b,df_d,m_gt=removedfault(df_b, df_d,m_gt)
    a=np.array(df_b["AD"])
    d=np.array(df_b["DP"])
    a=np.array(list(map(float,a)))
    d=np.array(list(map(float,d)))
    model =vireoSNP.VireoBulk(Ndonors)

    model.fit(a,d,m_gt,learn_theta=forceLearnGT_tag )
    return (model)
   
def preprocess_with_anno(bulk_data1,donor_data1,anno_data1,forceLearnGT_tag):
    bulk_data=load_VCF(bulk_data1, biallelic_only=True,load_sample=False,format_list=None,)
    donor_vcf=load_VCF( donor_data1,  biallelic_only=True,    load_sample=True, sparse=False , format_list=None)
    df_anno=pd.read_csv(anno_data1)
    Ndonors=len(donor_vcf['samples'])
    df_b,df_d,m_GT=rematch_merge_vcf_annotated(bulk_data,donor_vcf,df_anno)
    df_b,df_d,m_GT=removedfault(df_b,df_d,m_GT)
    a=np.array(df_b["AD"])
    d=np.array(df_b["DP"])
    a=np.array(list(map(float,a)))
    d=np.array(list(map(float,d)))
    model =vireoSNP.VireoBulk(Ndonors)
    gene_unique=np.unique(df_d["genes"].tolist())
    model.fit(a,d,m_GT,learn_theta=forceLearnGT_tag )
    model_return=model
    premodel=model.psi
    bulk_demuti_result= prediction_bulk(df_b,df_d,m_GT,donor_vcf,premodel,gene_unique)
    return(model_return,bulk_demuti_result)
    
main()
   
   
    
