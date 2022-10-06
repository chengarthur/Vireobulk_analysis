#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 20:25:39 2022

@author: zhizhiflyfly
"""

import argparse
from vireobulk import *
from version import __version__
import pandas as pd
import sys

def main():
   parser = argparse.ArgumentParser(description='input bulk pile up result reference,and annotation for SNPs(optional)')
   
   
   ##cell_data, the SNPs for bulk RNAseq data
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


   
   
   
   
    
if __name__ == "__main__":
    main()