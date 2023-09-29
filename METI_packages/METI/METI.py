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

def meta_gene_plot(img, 
            sudo_adata, 
            genes, 
            resize_factor, 
            binary, 
            res=50, 
            num_required=1, 
            pooling="min", 
            rseed=100, 
            tseed=100, 
            nseed=100, 
            nChannel=100, 
            threshold=1000, 
            radius=3, 
            minLabels=30, 
            train_refine=True, 
            plot_intermedium=False,
            target_size="small",
            min_UMI=5):
    #-------------------------------Image band-------------------------------------------------#
    if target_size=="small":
        target_ratio=2/3
    elif target_size=="large":
        target_ratio=1/3
    else:
        print("target_size not valid, Please specify (small / large).Use default small")
        target_ratio=1/2
    print("Computing image band...")
    resize_width=int(img.shape[1]*resize_factor)
    resize_height=int(img.shape[0]*resize_factor)
    img_resized=(img * np.dstack([(binary!=0)]*3)).astype(np.uint8)
    img_resized = cv2.resize(img_resized, (resize_width, resize_height), interpolation = cv2.INTER_AREA)
    gray_resized = cv2.cvtColor(img_resized,cv2.COLOR_BGR2GRAY)
    gray_resized=gray_resized.reshape(list(gray_resized.shape)+[1])
    binary_resized=cv2.resize(binary, (resize_width, resize_height), interpolation = cv2.INTER_AREA)
    img_band1=(gray_resized-np.min(gray_resized))/(np.max(gray_resized)-np.min(gray_resized))
    #-------------------------------Gene band-------------------------------------------------#
    print("Computing gene band...")
    genes=list(set([i for i in genes if i in sudo_adata.var.index ]))
    assert num_required<=len(genes)
    tmp=sudo_adata.X[:,sudo_adata.var.index.isin(genes)]
    tmp=(tmp-np.min(tmp, 0))/(np.max(tmp, 0)-np.min(tmp, 0))
    tmp = np.partition(tmp, -num_required, axis=1)[:,-num_required:] #Select top num_required
    if pooling=="mean":
        sudo_adata.obs["meta"]=np.mean(tmp, axis=1)
    elif pooling=="min":
        sudo_adata.obs["meta"]=np.min(tmp, axis=1)
    else:
        print("Error! Pooling logic not understood.")
    gene_img=np.zeros(list(img.shape[0:2])+[1])
    for _, row in sudo_adata.obs.iterrows():
        x=row["x"]
        y=row["y"]
        exp=row["meta"]
        gene_img[int(x-res/2):int(x+res/2), int(y-res/2):int(y+res/2),0]=exp
    gene_img[binary==0]=0
    return sudo_adata


def annotation(img, 
			sudo_adata, 
			genes, 
			resize_factor, 
			binary, 
			res=50, 
			num_required=1, 
			pooling="min", 
			rseed=100, 
			tseed=100, 
			nseed=100, 
			nChannel=100, 
			threshold=1000, 
			radius=3, 
			minLabels=30, 
			train_refine=True, 
			plot_intermedium=False,
			target_size="small",
			min_UMI=1):
	#-------------------------------Image band-------------------------------------------------#
	if target_size=="small":
		target_ratio=0.7
	elif target_size=="large":
		target_ratio=1/3
	else:
		print("target_size not valid, Please specify (small / large).Use default small")
		target_ratio=1/2
	print("Computing image band...")
	resize_width=int(img.shape[1]*resize_factor)
	resize_height=int(img.shape[0]*resize_factor)
	img_resized=(img * np.dstack([(binary!=0)]*3)).astype(np.uint8)
	img_resized = cv2.resize(img_resized, (resize_width, resize_height), interpolation = cv2.INTER_AREA)
	gray_resized = cv2.cvtColor(img_resized,cv2.COLOR_BGR2GRAY)
	gray_resized=gray_resized.reshape(list(gray_resized.shape)+[1])
	binary_resized=cv2.resize(binary, (resize_width, resize_height), interpolation = cv2.INTER_AREA)
	img_band1=(gray_resized-np.min(gray_resized))/(np.max(gray_resized)-np.min(gray_resized))
	#-------------------------------Gene band-------------------------------------------------#
	print("Computing gene band...")
	genes=list(set([i for i in genes if i in sudo_adata.var.index ]))
	assert num_required<=len(genes)
	tmp=sudo_adata.X[:,sudo_adata.var.index.isin(genes)]
	tmp=(tmp-np.min(tmp, 0))/(np.max(tmp, 0)-np.min(tmp, 0))
	tmp = np.partition(tmp, -num_required, axis=1)[:,-num_required:] #Select top num_required
	if pooling=="mean":
		sudo_adata.obs["meta"]=np.mean(tmp, axis=1)
	elif pooling=="min":
		sudo_adata.obs["meta"]=np.min(tmp, axis=1)
	else:
		print("Error! Pooling logic not understood.")
	gene_img=np.zeros(list(img.shape[0:2])+[1])
	for _, row in sudo_adata.obs.iterrows():
		x=row["x"]
		y=row["y"]
		exp=row["meta"]
		gene_img[int(x-res/2):int(x+res/2), int(y-res/2):int(y+res/2),0]=exp
	gene_img[binary==0]=0
	gene_img_scaled = (gene_img * 255).astype(np.uint8)
	gene_band1 = gene_img_scaled[:,:,0]
	gene_band1 = gene_band1/255
	#Filter on min UMI
# 	gene_band1[gene_band1<=np.log(min_UMI+1)]=0
	gene_band1=cv2.resize(gene_band1, (resize_width, resize_height), interpolation = cv2.INTER_AREA)
	gene_band1=(gene_band1-np.min(gene_band1))/(np.max(gene_band1)-np.min(gene_band1))
	gene_band1=gene_band1.reshape(list(gene_band1.shape)+[1])
	#-------------------------------TESLA-------------------------------------------------#
	print("Running TESLA...")
	assert np.max(img_band1)==1 and np.min(img_band1)==0 and np.max(gene_band1)==1 and np.min(gene_band1)==0
	data=np.concatenate((img_band1, gene_band1), axis=2)
	random.seed(rseed)
	torch.manual_seed(tseed)
	np.random.seed(nseed)
	from TESLA import TESLA
	from TESLA import refine_clusters
	tesla=TESLA()
	tesla.train(input=data, use_cuda=False, train_refine=train_refine, radius=radius, nChannel=nChannel, lr=0.1, 
		minLabels=minLabels, maxIter=30, stepsize_sim=1, stepsize_con=5, threshold=threshold, plot_intermedium=plot_intermedium, plot_dir=None)
	prob, pred = tesla.predict(data)
	pred_refined=refine_clusters(pred=pred, resize_height=resize_height, resize_width=resize_width, threshold=threshold, radius=radius)
	nLabels=len(np.unique(pred))
	mainLabels=len(np.unique(pred_refined))
	#-----------------------------------Find target cluster---------------------------#
	print("Finding target clusters...")	
	marker=gene_band1.flatten()
	clusters=pred_refined.flatten()
	c_m={} #cluster_marker expression
	for i in np.unique(clusters):
		c_m[i]=np.mean(marker[clusters==i])
	c_m = sorted(c_m.items(), key=lambda x: x[1], reverse=True)
	target_clusters=list(filter(lambda x: (x[1] > c_m[0][1]*target_ratio), c_m))
	target_clusters=[x[0] for x in target_clusters]
	print("c_m:\n", c_m, "\n", "Target clusters:\n", target_clusters, sep = '')
	return pred_refined, target_clusters, c_m





