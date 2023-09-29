
<h1><center>METI Tutorial</center></h1>

## Outline
1. [Installation]
2. [Import modules]
3. [Quality control]
4. [Read in data]

### 1. Installation
The installation should take a few minutes on a normal computer. 

Now you can install the current release of METI by the following three ways:
#### 1.1 PyPI: Directly install the package from PyPI.


```python
pip3 install METI
#Note: you need to make sure that the pip is for python3，or you can install METI by
python3 -m pip install METI==0.1
#If you do not have permission (when you get a permission denied error), you can install METI by
pip3 install --user METI==0.1
```
#### 1.2 Github
Download the package from Github and install it locally:


```python
git clone https://github.com/Flashiness/METI
cd METI/METI_package/
python3 setup.py install --user
```

### 2. Import modules

```python
import torch
import csv,re, time
import pickle
import random
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
from scipy import stats
from scipy.sparse import issparse
import scanpy as sc
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import cv2
import TESLA as tesla
from IPython.display import Image
import scipy.sparse
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scanpy import read_10x_h5
import PIL
from PIL import Image as IMAGE
import os
os.environ['KMP_DUPLICATE_LIB_OK']='True'
```

```python
METI.__version__
```
### 3. Quality control
```python
adata=sc.read_visium("/tutorial/data/Spaceranger/")
spatial=pd.read_csv("/tutorial/data/Spaceranger/tissue_positions_list.csv",sep=",",header=None,na_filter=False,index_col=0)

adata
```
    AnnData object with n_obs × n_vars = 3875 × 17943
    obs: 'in_tissue', 'array_row', 'array_col'
    var: 'gene_ids', 'feature_types', 'genome'
    uns: 'spatial'
    obsm: 'spatial'

```python
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

adata
```
    AnnData object with n_obs × n_vars = 3875 × 17943
    obs: 'in_tissue', 'array_row', 'array_col', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'pct_counts_in_top_50_genes', 'pct_counts_in_top_100_genes', 'pct_counts_in_top_200_genes', 'pct_counts_in_top_500_genes', 'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt'
    var: 'gene_ids', 'feature_types', 'genome', 'mt', 'n_cells_by_counts', 'mean_counts', 'log1p_mean_counts', 'pct_dropout_by_counts', 'total_counts', 'log1p_total_counts'
    uns: 'spatial'
    obsm: 'spatial'

```python
plt.rcParams["figure.figsize"] = (8, 8)
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"], size = 1.5, save = 'QC_plot.pdf')
```

**UMI counts**[embed](./sample_results/UMI_count.pdf)[/embed]


### 4. Read in data
The current version of METI requres three input data.
1. The gene expression matrix(n by k): expression_matrix.h5;
2. Spatial coordinateds of samplespositions.txt;
3. Histology image(optional): histology.tif, can be tif or png or jepg.

The gene expreesion data can be stored as an AnnData object. AnnData stores a data matrix .X together with annotations of observations .obs, variables .var and unstructured annotations .uns. 




















