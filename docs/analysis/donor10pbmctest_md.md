```python
## header file function (need to be script )
import sys
import numpy as np
import vireoSNP
from vireoSNP.utils.vireo_base import match
from vireoSNP.utils.vcf_utils import load_VCF, write_VCF, parse_donor_GPb
from vireoSNP.utils.vcf_utils import read_sparse_GeneINFO, GenoINFO_maker
import matplotlib.pyplot as plt
import pandas as pd
from vireoSNP.utils.io_utils import match_donor_VCF
import re as re
import os
```


```python
os.chdir("/home/ccheng/Vireobulk_analysis/")
```


```python
from base import merged_VCF_to_sdf,rematch_merge_vcf_annotated,find_all,removedfault,prediction_bulk
```


```python
##loadfile

donor_vcf=load_VCF(   "/media/wgs/subset_stimulation_bam/filter_pbmc10donors.vcf.gz",    biallelic_only=True,    load_sample=True, sparse=False , format_list=None)
####count based model
##bulk_merged_data=load_VCF(   "/media/wgs/subset_stimulation_bam/pbmc10donor_bulk/cellSNP.cells.vcf.gz",    biallelic_only=True,    load_sample=True,    format_list=None,)
#### UMI based model
bulk_merged_data=load_VCF("/media/wgs/subset_stimulation_bam/10pbmcdonor/cellSNP.base.vcf.gz", biallelic_only=True,load_sample=False,format_list=None,)
###processedcsv
df_anno=pd.read_csv("/media/wgs/subset_stimulation_bam/10pbmcdonor/pbmcvcf_annotated_processed.csv")

```


```python
#data clean
df_b,df_d,m_GT=rematch_merge_vcf_annotated(bulk_merged_data,donor_vcf,df_anno)
df_b,df_d,m_GT=removedfault(df_b,df_d,m_GT)
```


```python
df_b

```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>variants</th>
      <th>AD</th>
      <th>DP</th>
      <th>OTH</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1_14574_A_G</td>
      <td>2.0</td>
      <td>35.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1_14590_G_A</td>
      <td>0.0</td>
      <td>50.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1_14599_T_A</td>
      <td>0.0</td>
      <td>53.0</td>
      <td>2.0</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr1_14604_A_G</td>
      <td>1.0</td>
      <td>60.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>4</th>
      <td>chr1_14610_T_C</td>
      <td>1.0</td>
      <td>62.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>128180</th>
      <td>chrX_155999818_C_T</td>
      <td>66.0</td>
      <td>80.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>128181</th>
      <td>chrX_155999869_A_G</td>
      <td>88.0</td>
      <td>89.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>128182</th>
      <td>chrX_155999957_G_C</td>
      <td>73.0</td>
      <td>86.0</td>
      <td>0.0</td>
    </tr>
    <tr>
      <th>128183</th>
      <td>chrX_156025297_G_A</td>
      <td>65.0</td>
      <td>195.0</td>
      <td>1.0</td>
    </tr>
    <tr>
      <th>128187</th>
      <td>chrY_3120449_C_G</td>
      <td>0.0</td>
      <td>28.0</td>
      <td>1.0</td>
    </tr>
  </tbody>
</table>
<p>127155 rows × 4 columns</p>
</div>




```python
df_d
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>variants</th>
      <th>function</th>
      <th>genes</th>
      <th>avsnp</th>
      <th>GT</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>chr1_14574_A_G</td>
      <td>ncRNA_exonic</td>
      <td>WASH7P</td>
      <td>rs28503599</td>
      <td>[0/0, 0/1, 0/1, 0/0, 0/0, 0/0, 0/1, 0/0, 0/1, ...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>chr1_14590_G_A</td>
      <td>ncRNA_exonic</td>
      <td>WASH7P</td>
      <td>rs707679</td>
      <td>[0/0, 0/1, 0/1, 0/0, 0/0, 0/1, 0/1, 0/0, 0/0, ...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>chr1_14599_T_A</td>
      <td>ncRNA_exonic</td>
      <td>WASH7P</td>
      <td>rs707680</td>
      <td>[0/1, 0/1, 0/1, 0/0, 0/0, 0/1, 0/1, 0/1, 0/1, ...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>chr1_14604_A_G</td>
      <td>ncRNA_exonic</td>
      <td>WASH7P</td>
      <td>rs541940975</td>
      <td>[0/1, 0/1, 0/1, 0/0, 0/0, 0/1, 0/1, 0/1, 0/1, ...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>chr1_14610_T_C</td>
      <td>ncRNA_exonic</td>
      <td>WASH7P</td>
      <td>rs878986575</td>
      <td>[0/1, 0/1, 0/1, 0/0, 0/0, 0/1, 0/1, 0/1, 0/1, ...</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>128180</th>
      <td>chrX_155999818_C_T</td>
      <td>intronic</td>
      <td>IL9R</td>
      <td>rs3093470</td>
      <td>[1/1, 0/0, 1/1, 1/1, 0/1, 0/1, 0/1, 1/1, 1/1, ...</td>
    </tr>
    <tr>
      <th>128181</th>
      <td>chrX_155999869_A_G</td>
      <td>intronic</td>
      <td>IL9R</td>
      <td>rs3091261</td>
      <td>[1/1, 1/1, 1/1, 1/1, 1/1, 1/1, 1/1, 1/1, 1/1, ...</td>
    </tr>
    <tr>
      <th>128182</th>
      <td>chrX_155999957_G_C</td>
      <td>intronic</td>
      <td>IL9R</td>
      <td>rs3093472</td>
      <td>[1/1, 0/0, 1/1, 1/1, 0/1, 0/1, 0/1, 1/1, 1/1, ...</td>
    </tr>
    <tr>
      <th>128183</th>
      <td>chrX_156025297_G_A</td>
      <td>downstream</td>
      <td>DDX11L16</td>
      <td>rs185932868</td>
      <td>[0/1, 0/0, 0/1, 0/0, 0/1, 0/0, 0/0, 0/1, 0/1, ...</td>
    </tr>
    <tr>
      <th>128187</th>
      <td>chrY_3120449_C_G</td>
      <td>intergenic</td>
      <td>LINC00278;TGIF2LY</td>
      <td>rs778218266</td>
      <td>[0/1, 1/1, 0/1, 1/1, 0/0, 1/1, 1/1, 0/0, 0/1, ...</td>
    </tr>
  </tbody>
</table>
<p>127155 rows × 5 columns</p>
</div>




```python
gene_unique=np.unique(df_d["genes"].tolist())
```


```python
gene_unique=np.unique(df_d["genes"].tolist())
a=np.array(df_b["AD"])
d=np.array(df_b["DP"])
a=np.array(list(map(float,a)))
d=np.array(list(map(float,d)))
```


```python
model =vireoSNP.VireoBulk(10)
```


```python
model.fit(a,d,m_GT,learn_theta=True )
```


```python
model.theta
```




    array([0.01164797, 0.44852633, 0.99267153])




```python
donor_vcf['samples']
```




    ['18-G-017_R1',
     '18-G-018_R1',
     '18-G-019_R1',
     '18-G-020_R1',
     '18-G-021_R1',
     '18-G-022_R1',
     '18-G-023_R1',
     '18-G-024_R1',
     '18-G-025_R1',
     '18-G-026_R1']




```python
model.psi
```




    array([0.06418119, 0.04299516, 0.05283177, 0.08967812, 0.14387843,
           0.06906367, 0.04872266, 0.11876265, 0.29342823, 0.07645813])




```python
premodel=model.psi
```


```python
###gene model

bulk_demuti_result= prediction_bulk(df_b,df_d,m_GT,donor_vcf,premodel,gene_unique)
```


```python
###the 10 donor gene level demutiplexing may not as accurate as it in two donors  this result is just for demonstration
bulk_demuti_result
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>18-G-017_R1</th>
      <th>18-G-018_R1</th>
      <th>18-G-019_R1</th>
      <th>18-G-020_R1</th>
      <th>18-G-021_R1</th>
      <th>18-G-022_R1</th>
      <th>18-G-023_R1</th>
      <th>18-G-024_R1</th>
      <th>18-G-025_R1</th>
      <th>18-G-026_R1</th>
      <th>chi</th>
      <th>p</th>
      <th>numSNP</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>A1BG</th>
      <td>0.240487</td>
      <td>0.031524</td>
      <td>0.246912</td>
      <td>0.0424</td>
      <td>0.212023</td>
      <td>0.01513</td>
      <td>0.073156</td>
      <td>0.081881</td>
      <td>0.022049</td>
      <td>0.034438</td>
      <td>3.290594</td>
      <td>0.951652</td>
      <td>3</td>
    </tr>
    <tr>
      <th>A1BG-AS1</th>
      <td>0.322538</td>
      <td>0.000538</td>
      <td>0.043626</td>
      <td>0.000298</td>
      <td>0.024924</td>
      <td>0.110349</td>
      <td>0.003099</td>
      <td>0.307019</td>
      <td>0.003918</td>
      <td>0.183691</td>
      <td>2.721691</td>
      <td>0.974315</td>
      <td>2</td>
    </tr>
    <tr>
      <th>A2M</th>
      <td>0.111379</td>
      <td>0.186228</td>
      <td>0.015524</td>
      <td>0.070369</td>
      <td>0.002524</td>
      <td>0.00208</td>
      <td>0.113103</td>
      <td>0.210372</td>
      <td>0.260654</td>
      <td>0.027765</td>
      <td>6.853791</td>
      <td>0.652339</td>
      <td>2</td>
    </tr>
    <tr>
      <th>A2M;PZP</th>
      <td>0.297699</td>
      <td>0.057649</td>
      <td>0.159085</td>
      <td>0.049461</td>
      <td>0.025333</td>
      <td>0.013277</td>
      <td>0.182469</td>
      <td>0.030398</td>
      <td>0.16116</td>
      <td>0.023469</td>
      <td>5.948257</td>
      <td>0.745086</td>
      <td>1</td>
    </tr>
    <tr>
      <th>A2MP1</th>
      <td>0.054614</td>
      <td>0.433896</td>
      <td>0.016346</td>
      <td>0.109016</td>
      <td>0.105064</td>
      <td>0.050424</td>
      <td>0.014692</td>
      <td>0.170016</td>
      <td>0.045928</td>
      <td>0.000005</td>
      <td>22.90447</td>
      <td>0.006414</td>
      <td>3</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>ZXDC</th>
      <td>0.183335</td>
      <td>0.069798</td>
      <td>0.073204</td>
      <td>0.156388</td>
      <td>0.011695</td>
      <td>0.013849</td>
      <td>0.084582</td>
      <td>0.1278</td>
      <td>0.164845</td>
      <td>0.114505</td>
      <td>2.391586</td>
      <td>0.983658</td>
      <td>13</td>
    </tr>
    <tr>
      <th>ZYG11B</th>
      <td>0.084576</td>
      <td>0.162609</td>
      <td>0.022348</td>
      <td>0.183719</td>
      <td>0.0511</td>
      <td>0.027481</td>
      <td>0.088382</td>
      <td>0.223326</td>
      <td>0.086492</td>
      <td>0.069967</td>
      <td>10.147248</td>
      <td>0.338697</td>
      <td>9</td>
    </tr>
    <tr>
      <th>ZYX</th>
      <td>0.026771</td>
      <td>0.061399</td>
      <td>0.018724</td>
      <td>0.008064</td>
      <td>0.084577</td>
      <td>0.122853</td>
      <td>0.104039</td>
      <td>0.288391</td>
      <td>0.283356</td>
      <td>0.001826</td>
      <td>7.684969</td>
      <td>0.566179</td>
      <td>2</td>
    </tr>
    <tr>
      <th>ZZEF1</th>
      <td>0.063709</td>
      <td>0.064758</td>
      <td>0.008241</td>
      <td>0.103454</td>
      <td>0.050555</td>
      <td>0.040779</td>
      <td>0.000061</td>
      <td>0.09356</td>
      <td>0.474027</td>
      <td>0.100855</td>
      <td>1.673724</td>
      <td>0.995641</td>
      <td>29</td>
    </tr>
    <tr>
      <th>ZZZ3</th>
      <td>0.036555</td>
      <td>0.076785</td>
      <td>0.262047</td>
      <td>0.230531</td>
      <td>0.154715</td>
      <td>0.06554</td>
      <td>0.002896</td>
      <td>0.050482</td>
      <td>0.036244</td>
      <td>0.084205</td>
      <td>4.397809</td>
      <td>0.883336</td>
      <td>25</td>
    </tr>
  </tbody>
</table>
<p>13350 rows × 13 columns</p>
</div>




```python

```
