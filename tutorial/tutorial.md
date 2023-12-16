
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
import METI as meti
import tifffile
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
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"], size = 1.5, save = 'UMI_count.png')
```

**UMI counts**![](./sample_results/UMI_count_Scanpy_54078.png)


### 4. Read in data
The current version of METI requres three input data.
1. The gene expression matrix(n by k): expression_matrix.h5;
2. Spatial coordinateds of samplespositions.txt;
3. Histology image(optional): histology.tif, can be tif or png or jepg.

```python
#================== 3. Read in data ==================#
#Read original 10x_h5 data and save it to h5ad
from scanpy import read_10x_h5
adata = read_10x_h5("../tutorial/data/filtered_feature_bc_matrix.h5")
spatial=pd.read_csv("../tutorial/data/tissue_positions_list.csv",sep=",",header=None,na_filter=False,index_col=0) 
adata.obs["x1"]=spatial[1]
adata.obs["x2"]=spatial[2]
adata.obs["x3"]=spatial[3]
adata.obs["x4"]=spatial[4]
adata.obs["x5"]=spatial[5] 
adata.obs["array_x"]=adata.obs["x2"]
adata.obs["array_y"]=adata.obs["x3"]
adata.obs["pixel_x"]=adata.obs["x4"]
adata.obs["pixel_y"]=adata.obs["x5"]
#Select captured samples
adata=adata[adata.obs["x1"]==1]
adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")
adata.write_h5ad("../tutorial/data/1957495_data.h5ad")

#Read in gene expression and spatial location
counts=sc.read("../tutorial/data/1957495_data.h5ad")
#Read in hitology image
PIL.Image.MAX_IMAGE_PIXELS = None
img = IMAGE.open(r"../tutorial/data/histology.tif") 
img = np.array(img)

#if your image has 4 dimensions, only keep first 3 dims
img=img[...,:3]
```

### 5. Gene expression enhancement
#### 5.1 Preprocessing
```python
resize_factor=1000/np.min(img.shape[0:2])
resize_width=int(img.shape[1]*resize_factor)
resize_height=int(img.shape[0]*resize_factor)
counts.var.index=[i.upper() for i in counts.var.index]
counts.var_names_make_unique()
counts.raw=counts
sc.pp.log1p(counts) # impute on log scale
if issparse(counts.X):counts.X=counts.X.A
```
#### 5.2 Contour detection
```python
# Detect contour using cv2
cnt=tesla.cv2_detect_contour(img, apertureSize=5,L2gradient = True)

binary=np.zeros((img.shape[0:2]), dtype=np.uint8)
cv2.drawContours(binary, [cnt], -1, (1), thickness=-1)
#Enlarged filter
cnt_enlarged = tesla.scale_contour(cnt, 1.05)
binary_enlarged = np.zeros(img.shape[0:2])
cv2.drawContours(binary_enlarged, [cnt_enlarged], -1, (1), thickness=-1)
img_new = img.copy()
cv2.drawContours(img_new, [cnt], -1, (255), thickness=20)
img_new=cv2.resize(img_new, ((resize_width, resize_height)))
cv2.imwrite('../tutorial/data/cnt_1957495.jpg', img_new)
Image(filename='../tutorial/data/cnt_1957495.jpg')
```

#### 5.3 Gene expression enhancement
```python
#Set size of superpixel
res=40
# Note, if the numer of superpixels is too large and take too long, you can increase the res to 100
enhanced_exp_adata=tesla.imputation(img=img, raw=counts, cnt=cnt, genes=counts.var.index.tolist(), shape="None", res=res, s=1, k=2, num_nbs=10)
enhanced_exp_adata.write_h5ad("../tutorial/data/enhanced_exp.h5ad")
```
    Total number of sudo points:  92844
    Calculating spot 0
    Calculating spot 1000
    Calculating spot 2000
    Calculating spot 3000
    Calculating spot 4000
    Calculating spot 5000
    Calculating spot 6000
    Calculating spot 7000
    Calculating spot 8000
    Calculating spot 9000
    Calculating spot 10000
    Calculating spot 11000
    Calculating spot 12000
    Calculating spot 13000
    Calculating spot 14000
    Calculating spot 15000
    Calculating spot 16000
    Calculating spot 17000
    Calculating spot 18000
    Calculating spot 19000
    Calculating spot 20000
    Calculating spot 21000
    Calculating spot 22000
    Calculating spot 23000
    Calculating spot 24000
    Calculating spot 25000
    Calculating spot 26000
    Calculating spot 27000
    Calculating spot 28000
    Calculating spot 29000
    Calculating spot 30000
    Calculating spot 31000
    Calculating spot 32000
    Calculating spot 33000
    Calculating spot 34000
    Calculating spot 35000
    Calculating spot 36000
    Calculating spot 37000
    Calculating spot 38000
    Calculating spot 39000
    Calculating spot 40000
    Calculating spot 41000
    Calculating spot 42000
    Calculating spot 43000
    Calculating spot 44000
    Calculating spot 45000
    Calculating spot 46000
    Calculating spot 47000
    Calculating spot 48000
    Calculating spot 49000
    Calculating spot 50000
    Calculating spot 51000
    Calculating spot 52000
    Calculating spot 53000
    Calculating spot 54000
    Calculating spot 55000
    Calculating spot 56000
    Calculating spot 57000
    Calculating spot 58000
    Calculating spot 59000
    Calculating spot 60000
    Calculating spot 61000
    Calculating spot 62000
    Calculating spot 63000
    Calculating spot 64000
    Calculating spot 65000
    Calculating spot 66000
    Calculating spot 67000
    Calculating spot 68000
    Calculating spot 69000
    Calculating spot 70000
    Calculating spot 71000
    Calculating spot 72000
    Calculating spot 73000
    Calculating spot 74000
    Calculating spot 75000
    Calculating spot 76000
    Calculating spot 77000
    Calculating spot 78000
    Calculating spot 79000
    Calculating spot 80000
    Calculating spot 81000
    Calculating spot 82000
    Calculating spot 83000
    Calculating spot 84000
    Calculating spot 85000
    Calculating spot 86000
    Calculating spot 87000
    Calculating spot 88000
    Calculating spot 89000
    Calculating spot 90000
    Calculating spot 91000
    Calculating spot 92000
    --- 280.4760136604309 seconds ---
    Imputing spot 0
    Imputing spot 1000
    Imputing spot 2000
    Imputing spot 3000
    Imputing spot 4000
    Imputing spot 5000
    Imputing spot 6000
    Imputing spot 7000
    Imputing spot 8000
    Imputing spot 9000
    Imputing spot 10000
    Imputing spot 11000
    Imputing spot 12000
    Imputing spot 13000
    Imputing spot 14000
    Imputing spot 15000
    Imputing spot 16000
    Imputing spot 17000
    Imputing spot 18000
    Imputing spot 19000
    Imputing spot 20000
    Imputing spot 21000
    Imputing spot 22000
    Imputing spot 23000
    Imputing spot 24000
    Imputing spot 25000
    Imputing spot 26000
    Imputing spot 27000
    Imputing spot 28000
    Imputing spot 29000
    Imputing spot 30000
    Imputing spot 31000
    Imputing spot 32000
    Imputing spot 33000
    Imputing spot 34000
    Imputing spot 35000
    Imputing spot 36000
    Imputing spot 37000
    Imputing spot 38000
    Imputing spot 39000
    Imputing spot 40000
    Imputing spot 41000
    Imputing spot 42000
    Imputing spot 43000
    Imputing spot 44000
    Imputing spot 45000
    Imputing spot 46000
    Imputing spot 47000
    Imputing spot 48000
    Imputing spot 49000
    Imputing spot 50000
    Imputing spot 51000
    Imputing spot 52000
    Imputing spot 53000
    Imputing spot 54000
    Imputing spot 55000
    Imputing spot 56000
    Imputing spot 57000
    Imputing spot 58000
    Imputing spot 59000
    Imputing spot 60000
    Imputing spot 61000
    Imputing spot 62000
    Imputing spot 63000
    Imputing spot 64000
    Imputing spot 65000
    Imputing spot 66000
    Imputing spot 67000
    Imputing spot 68000
    Imputing spot 69000
    Imputing spot 70000
    Imputing spot 71000
    Imputing spot 72000
    Imputing spot 73000
    Imputing spot 74000
    Imputing spot 75000
    Imputing spot 76000
    Imputing spot 77000
    Imputing spot 78000
    Imputing spot 79000
    Imputing spot 80000
    Imputing spot 81000
    Imputing spot 82000
    Imputing spot 83000
    Imputing spot 84000
    Imputing spot 85000
    Imputing spot 86000
    Imputing spot 87000
    Imputing spot 88000
    Imputing spot 89000
    Imputing spot 90000
    Imputing spot 91000
    Imputing spot 92000

### 6. Goblet marker gene expression
```python
#================ determine if markers are in ===============#
enhanced_exp_adata=sc.read("..tutorial/data/enhanced_exp.h5ad")
markers = ["MS4A10", "MGAM", "CYP4F2", "XPNPEP2", "SLC5A9", "SLC13A2", "SLC28A1", "MEP1A", "ABCG2", "ACE2"]
for i in range(len(markers)):
    if markers[i] in enhanced_exp_adata.var.index: print("yes")
    else: print(markers[i])
```

```python
save_dir="..tutorial/data/Goblet/"
if not os.path.exists(save_dir):os.mkdir(save_dir)
```

```python
#================ Plot gene expression image ===============#
markers = ["MS4A10", "MGAM", "CYP4F2", "XPNPEP2", "SLC5A9", "SLC13A2", "SLC28A1", "MEP1A", "ABCG2", "ACE2"]
for i in range(len(markers)):
    cnt_color = clr.LinearSegmentedColormap.from_list('magma', ["#000003",  "#3b0f6f",  "#8c2980",   "#f66e5b", "#fd9f6c", "#fbfcbf"], N=256)
    g=markers[i]
    enhanced_exp_adata.obs[g]=enhanced_exp_adata.X[:,enhanced_exp_adata.var.index==g]
    fig=sc.pl.scatter(enhanced_exp_adata,alpha=1,x="y",y="x",color=g,color_map=cnt_color,show=False,size=10)
    fig.set_aspect('equal', 'box')
    fig.invert_yaxis()
    plt.gcf().set_dpi(600)
    fig.figure.show()

    plt.savefig(save_dir + str(markers[i]) + ".png", dpi=600)
    plt.close()
```
**Gene expression**![](./sample_results/ABCG2.png)

```python
#================ Plot meta gene expression image ===============#
enhanced_exp_adata=sc.read("/Users/jjiang6/Desktop/UTH/MDA GRA/Spatial transcriptome/Cell Segmentation/With Jian Hu/S1_54078/TESLA/enhanced_exp.h5ad")
genes =   ["MS4A10", "MGAM", "CYP4F2", "XPNPEP2", "SLC5A9", "SLC13A2", "SLC28A1", "MEP1A", "ABCG2", "ACE2"]
    
sudo_adata = meti.meta_gene_plot(img=img, 
                                binary=binary,
                                sudo_adata=enhanced_exp_adata, 
                                genes=genes, 
                                resize_factor=resize_factor,
                                target_size="small")

cnt_color = clr.LinearSegmentedColormap.from_list('magma', ["#000003",  "#3b0f6f",  "#8c2980",   "#f66e5b", "#fd9f6c", "#fbfcbf"], N=256)
fig=sc.pl.scatter(sudo_adata,alpha=1,x="y",y="x",color='meta',color_map=cnt_color,show=False,size=5)
fig.set_aspect('equal', 'box')
fig.invert_yaxis()
plt.gcf().set_dpi(600)
fig.figure.show()

plt.savefig(save_dir + "Goblet_meta.png", dpi=600)
plt.close()
```
**Meta gene expression**![](./sample_results/Goblet_meta.png)

### 7. Region annotation
```python
genes=["MS4A10", "MGAM", "CYP4F2", "XPNPEP2", "SLC5A9", "SLC13A2", "SLC28A1", "MEP1A", "ABCG2", "ACE2"]
genes=list(set([i for i in genes if i in enhanced_exp_adata.var.index ]))
#target_size can be set to "small" or "large".
pred_refined, target_clusters, c_m=meti.annotation(img=img, 
                                                    binary=binary,
                                                    sudo_adata=enhanced_exp_adata, 
                                                    genes=genes, 
                                                    resize_factor=resize_factor,
                                                    num_required=1, 
                                                    target_size="small")
#Plot
ret_img=tesla.visualize_annotation(img=img, 
                              binary=binary, 
                              resize_factor=resize_factor,
                              pred_refined=pred_refined, 
                              target_clusters=target_clusters, 
                              c_m=c_m)

cv2.imwrite(save_dir + 'IME.jpg', ret_img)
Image(filename=save_dir + 'IME.jpg')
```
**Goblet region annotation**![](./sample_results/Goblet.jpg)


### 8. Segmentation
```python
img = tifffile.imread(r"../tutorial/data/1415785-6 Bx2.tif") 

patch_size = 4000
d0=int(np.ceil(img.shape[0]/patch_size)*patch_size)
d1=int(np.ceil(img.shape[1]/patch_size)*patch_size)

img_extended=np.concatenate((img, np.zeros((img.shape[0], d1-img.shape[1], 3), dtype=np.uint16)), axis=1)
img_extended=np.concatenate((img_extended, np.zeros((d0-img.shape[0], d1, 3), dtype=np.uint16)), axis=0)

#=====================================Split into patched=====================================================
x=[i*patch_size for i in range(int(img_extended.shape[0]/patch_size))]
x=np.repeat(x,int(img_extended.shape[1]/patch_size))
y=[i*patch_size for i in range(int(img_extended.shape[1]/patch_size))]
y=np.array(y*int(img_extended.shape[0]/patch_size))
patches=np.zeros((len(x), patch_size, patch_size, 3), dtype=np.uint16) #n*patch_size*patch_size*3

patch_info_list = [{'x': x, 'y': y} for x, y in zip(x, y)]
patch_info = pd.DataFrame(patch_info_list)
patch_info.to_csv('../tutorial/data/seg_results/patch_info.csv', index=False)

counter=0
for i in range(len(x)):
	x_tmp=int(x[i])
	y_tmp=int(y[i])
	patches[counter, :, :, :]=img_extended[x_tmp:x_tmp+patch_size,y_tmp:y_tmp+patch_size, :]
	counter+=1

#patch_info = pd.read_csv('../tutorial/data/seg_results/patch_info.csv')
save_dir = '../tutorial/data/seg_results/'

#=================================Kmeans Segmentation===================================
meti.Segment_Patches(patches, save_dir, n_clusters=10)
```
    Doing:  0 / 9
    Doing:  1 / 9
    Doing:  2 / 9
    Doing:  3 / 9
    Doing:  4 / 9
    Doing:  5 / 9
    Doing:  6 / 9
    Doing:  7 / 9
    Doing:  8 / 9

```python
pred_file_locs=[save_dir+"/patch"+str(j)+"_pred.npy" for j in range(patch_info.shape[0])]
dic_list=meti.get_color_dic(patches, seg_dir=save_dir)
masks_index=meti.Match_Masks(dic_list, num_mask_each=5, mapping_threshold1=30, mapping_threshold2=60)
```
    Doing  0
    Adding,  no match counts: 0 new rate: 1.0
    Adding,  no match counts: 1 new rate: 1.0
    Not adding,  no match counts: 4 new rate: 1.0
    Doing  1
    Not adding,  no match counts: 0 new rate: 0.6842105263157895
    Adding,  no match counts: 0 new rate: 1.0
    Not adding,  no match counts: 2 new rate: 0.6
    Doing  2
    Doing  3
    Doing  4
    Not adding,  no match counts: 0 new rate: 0.12903225806451613
    Adding,  no match counts: 1 new rate: 0.9
    Doing  5
    Doing  6
    Doing  7
    Doing  8
    
```python
masks=meti.Extract_Masks(masks_index, pred_file_locs, patch_size)
```
    Extracting mask  0
    Extracting mask  1
    Extracting mask  2
    Extracting mask  3

```python
combined_masks=meti.Combine_Masks(masks, patch_info, img.shape[0], img.shape[1])
```
    Combining mask  0
    Combining mask  1
    Combining mask  2
    Combining mask  3

```python
plot_dir = '../tutorial/data/seg_results/mask'

for i in range(masks.shape[0]): #Each mask
	print("Plotting mask ", str(i))
	ret=(combined_masks[i]*255)
	np.save(plot_dir + '/masks' + str(i) + '.npy', ret)
	cv2.imwrite(plot_dir+'/mask'+str(i)+'.png', ret.astype(np.uint8))
```
    Plotting mask  0
    Plotting mask  1
    Plotting mask  2
    Plotting mask  3

```python
#===========================find a correct mask meeting with your requirement, ex: nucleis channel
mask_nuclei = np.load('../tutorial/data/seg_results/masks4_nuclei.npy')
mask_nuclei = mask_nuclei.astype(np.uint8)
ret, labels = cv2.connectedComponents(mask_nuclei)
cc_features=Extract_CC_Features_each_CC(labels)

filtered_cc_index=cc_features[(cc_features["cc_areas"]>40) &(cc_features["area_ratios"]>0.5) & (cc_features["hw_ratios"]<3)]["cc_index"].tolist()

tmp=ret*(np.isin(ret, filtered_cc_index))
tmp=plot_cc(tmp)
cv2.imwrite(plot_dir+'/nuclei_filtered.png', tmp)
```

**nuclei segmentation**![](./sample_results/nuclei_filtered1.png)

**goblet filter**![](./sample_results/nuclei_filtered2.png)











