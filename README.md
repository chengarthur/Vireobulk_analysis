# Vireobulk： Bayesian methods for bulk multiplexed RNA demultiplexing 



Vireobulk is a tool that leverages SNPs to devolve mixed donor RNAseq samples, which includes two main  modules :

## [**Documents**](https://vireobulk-analysis.readthedocs.io/)

## Donor decomposition 
With given  donor SNP reference (VCF file) and piled up Bulk RNAseq result(VCF file ), Vireobulk can accurate estimates the percentage of in an mixed donor sample in situation such as bone marrow transplantation, chimeric organoids, pooled iPSC cell lines. In most cases, vireo bulk can demultiplex  
from 2 to more than 20 donors 
## Gene abundance inference

with detailed Gene-SNP annotation,  vireo bulk can even give the estimation of gene abundance of different donor 
**reminder**: the power of vireo-bulk rely on sequencing depth and the number of donor-Genes pairs.

## quick start 
### installation 
```


git clone https://github.com/chengarthur/Vireobulk_analysis.git
cd Vireobulk_analysis/
pip install -r requirements.txt
```

### run test1 Donor decomposition
```

python vireobulk.py -c "data/cellSNP.base.vcf.gz" -n test1 -d "data/filter_pbmc10donors.vcf.gz"
````
This step will generate summary file for demultiplexed donor ratio and model.theta file
demultiplexed donor ratio contains sample X ratio matrix in piled up sample 
mode theta represents the inferred B Allele frequency of model, the more closed to [0,0.5,1]，the more reliable of the result, usually BAF > 0.40 represent confident inference 


### run test1 gene Decomposition 
```
python vireobulk.py -c "data/cellSNP.base.vcf.gz" -n test1 - -n test1 -d "data/filter_pbmc10donors.vcf.gz" -a "data/scvcf_annotated_processed.csv" --gene



```
This step will take few minutes generate summary file for demultiplexed donor ratio and model.theta file ,and a gene X donor ratio file ,and withe the pvalue of if the Gene is differentially expressed among donors.





## support
The document  and detailed analysis, data preparation usage can be checked on documents[read the doc ] if you have any questions and issues to ask , please leave issue on GitHub or email  author @Cheng Chen  <zhizff74@connect.hku.hk>



