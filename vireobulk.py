#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 19:46:07 2022

@author: zhizhiflyfly
"""

## header file function (need to be script )
import sys
import numpy as np
import vireoSNP
from vireoSNP.utils.vcf_utils import read_sparse_GeneINFO, GenoINFO_maker
import pandas as pd

import re as re
def VCF_to_sdf(bulk_data):
    df_empty=pd.DataFrame(columns=["variants","AD","DP","GT"])
    bulk_dat = read_sparse_GeneINFO(bulk_data['GenoINFO'], keys=['AD', 'DP'])
    df_empty["variants"]=bulk_data["variants"]
    df_empty["AD"]=bulk_dat["AD"].toarray()
    df_empty["DP"]=bulk_dat["DP"].toarray()
    return (df_empty)

def rematch_vcf(bulk_data,donor_vcf):
    df_bulk=VCF_to_sdf(bulk_data)
    df_donor=pd.DataFrame(columns=["variants","GT"])
    df_donor["variants"]=donor_vcf["variants"]
    df_donor["GT"]=donor_vcf["GenoINFO"]["GT"]     
    donor_tensor = vireoSNP.vcf.parse_donor_GPb(donor_vcf['GenoINFO']["GT"], "GT")
    itersect_variants=pd.Series(list(set(df_donor["variants"]).intersection(set(df_bulk["variants"]))))
    id_d=df_donor["variants"].isin(itersect_variants)
    id_b=df_bulk["variants"].isin(itersect_variants)
    df_bulk=df_bulk[id_b]
    df_donor=df_donor[id_d]
    donor_tensor=donor_tensor[id_d]
    return(df_bulk,df_donor,donor_tensor)
def rematch_vcf_annotated(bulk_data,donor_vcf,donor_vcf_anno):
    df_bulk=VCF_to_sdf(bulk_data)
    df_donora=pd.DataFrame(columns=["variants","function","genes","avsnp","GT"])
    df_donora["variants"]=donor_vcf_anno["variants"]  
    df_donora["function"]=donor_vcf_anno["function."]
    df_donora["genes"]=donor_vcf_anno["genes"]
    df_donora["avsnp"]=donor_vcf_anno["avsnp"]
    df_donor=pd.DataFrame(columns=["variants","GT"])
    df_donor["variants"]=donor_vcf["variants"]
    df_donor["GT"]=donor_vcf["GenoINFO"]["GT"]     
    donor_tensor = vireoSNP.vcf.parse_donor_GPb(donor_vcf['GenoINFO']["GT"], "GT")
    itersect_variants=pd.Series(list(set(df_donora["variants"]).intersection(set(df_bulk["variants"]))))
    id_d=df_donor["variants"].isin(itersect_variants)
    id_b=df_bulk["variants"].isin(itersect_variants)
    id_da=df_donora["variants"].isin(itersect_variants)
    df_bulk=df_bulk[id_b]
    df_donor=df_donor[id_d]
    df_donora=df_donora[id_da]
    donor_tensor=donor_tensor[id_d]
    df_donora["GT"]=df_donora["GT"]=df_donor["GT"].tolist()
    return(df_bulk,df_donora,donor_tensor)
def merged_VCF_to_sdf(bulk_data_merge):
    df_empty=pd.DataFrame(columns=["variants","AD","DP","OTH"])
    df_empty["variants"]=bulk_data_merge["variants"]
    c=np.array(bulk_data_merge['FixedINFO']["INFO"])
    for i in range(len(c)):
        tem=c[i].tolist().split(";",-1)
        df_empty["AD"][i]=float(re.sub("AD=","",tem[0]))
        df_empty["DP"][i]=float(re.sub("DP=","",tem[1]))
        df_empty["OTH"][i]=float(re.sub("OTH=","",tem[2]))

    return(df_empty)
def rematch_merge_vcf_annotated(bulk_data,donor_vcf,donor_vcf_anno):
    df_bulk=merged_VCF_to_sdf(bulk_data)
    df_donora=pd.DataFrame(columns=["variants","function","genes","avsnp","GT"])
    df_donora["variants"]=donor_vcf_anno["variants"]  
    df_donora["function"]=donor_vcf_anno["function."]
    df_donora["genes"]=donor_vcf_anno["genes"]
    df_donora["avsnp"]=donor_vcf_anno["avsnp"]
    df_donor=pd.DataFrame(columns=["variants","GT"])
    df_donor["variants"]=donor_vcf["variants"]
    df_donor["GT"]=donor_vcf["GenoINFO"]["GT"]     
    donor_tensor = vireoSNP.vcf.parse_donor_GPb(donor_vcf['GenoINFO']["GT"], "GT")
    itersect_variants=pd.Series(list(set(df_donora["variants"]).intersection(set(df_bulk["variants"]))))
    id_d=df_donor["variants"].isin(itersect_variants)
    id_b=df_bulk["variants"].isin(itersect_variants)
    id_da=df_donora["variants"].isin(itersect_variants)
    df_bulk=df_bulk[id_b]
    df_donor=df_donor[id_d]
    df_donora=df_donora[id_da]
    donor_tensor=donor_tensor[id_d]
    df_donora["GT"]=df_donora["GT"]=df_donor["GT"].tolist()
    return(df_bulk,df_donora,donor_tensor)
def find_all(data, s):
    r_list = []
    for r in range(len(data)):
        r_list.append( s in data[r])
 
    return(r_list)
def removedfault(df_bulk,df_donora,donor_tensor):
    id_d=[]
    for bool in find_all(np.array(df_donora["GT"]).tolist(),"./."):
        id_d.append(not bool)
    df_donor=df_donora[id_d]
    df_gt=donor_tensor[id_d]
    id_b=df_bulk["variants"].isin(df_donor["variants"])
    df_b=df_bulk[id_b]
    return(df_b,df_donor,df_gt)

def prediction_bulk(df_b,df_d,df_gt,donor_vcf,pre_model,genelist):
    gene_unique=genelist
    new=pd.DataFrame(index=gene_unique,columns=donor_vcf['samples']+["chi","p","numSNP"])
    pre_psi=pre_model
    test_d=df_d
    test_b=df_b
    test_gt=df_gt
    sample_list=donor_vcf['samples']
    nk=len(donor_vcf['samples'])
    for i in range(len(gene_unique)):
        sub_d=test_d[test_d["genes"]==gene_unique[i]]
        sub_b=test_b[test_b["variants"].isin(test_d[test_d["genes"]==gene_unique[i]]["variants"])]
        sub_gt=test_gt[test_d["genes"]==gene_unique[i]]
        a=np.array(sub_b["AD"])
        d=np.array(sub_b["DP"])
        a=np.array(list(map(float,a)))
        d=np.array(list(map(float,d)))
        model =vireoSNP.VireoBulk(nk)
        model.fit(a,d,sub_gt,learn_theta=False )
        pv=vireoSNP.LikRatio_test(model.psi,pre_psi,a,d,sub_gt,model.theta)
        for j in range(len(model.psi)):
            new[sample_list[j]][i]=model.psi[j]
        new["chi"][i]=pv[0]
        new["p"][i]=pv[1]
        new["numSNP"][i]=len(a)
    return(new)