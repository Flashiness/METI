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

def patch_split_for_ST(img, patch_size, spot_info, x_name="pixel_x", y_name="pixel_y"):
	assert patch_size%2==0
	patches=np.zeros((spot_info.shape[0], patch_size, patch_size, 3), dtype=np.uint8) 
	counter=0
	for _, row in spot_info.iterrows():
		x_tmp=int(row[x_name])
		y_tmp=int(row[y_name])
		patches[counter, :, :, :]=img[int(x_tmp-patch_size/2):int(x_tmp+patch_size/2),int(y_tmp-patch_size/2):int(y_tmp+patch_size/2), :]
		filename = os.path.join('/rsrch4/home/genomic_med/jjiang6/Project1/S1_54078/Segmentation/NC_review_Goblet_seg/patches/', f"patch_{counter}.png")
		cv2.imwrite(filename, cv2.cvtColor(patches[counter, :, :, :], cv2.COLOR_RGB2BGR))
		counter+=1
	return patches

def Segment_Patches(patches, save_dir="./seg_results",n_clusters=10):
	patch_size=patches.shape[1]
	for i in range(patches.shape[0]):
		print("Doing: ", i, "/", patches.shape[0])
		patch=patches[i]
		pixel_values = patch.reshape((-1, 3))
		#convert to float
		pixel_values = np.float32(pixel_values)
		seed=100
		random.seed(seed)
		np.random.seed(seed)
		kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(pixel_values)
		pred=kmeans.labels_
		#centroids=kmeans.cluster_centers_
		np.save(save_dir+"/patch"+str(i)+"_pred.npy", pred)

def get_color_dic(patches, seg_dir,  refine=True):
	patch_index=list(range(patches.shape[0]))
	patch_size=patches.shape[1]
	dic_list=[]
	for i in patch_index:
		pred_refined=np.load(seg_dir+"/patch"+str(i)+"_pred.npy")
		patch=patches[i].copy()
		c_m={} #cluster_marker expression
		c_a=[(j, np.sum(pred_refined==j)) for j in np.unique(pred_refined)] #cluster_area
		c_a=sorted(c_a, key=lambda x: x[1], reverse=True)
		c_a=dict((x, y) for x, y in c_a)
		clusters=pred_refined.reshape(patch_size, patch_size)
		for j,  area_ratio in c_a.items():
			c_m[j]=np.median(patch[clusters==j], 0).tolist()
		dic_list.append(c_m)
	return dic_list

def Match_Masks(dic_list, num_mask_each=10, mapping_threshold1=30, mapping_threshold2=60, min_no_mach_rate=0.2, min_unused_channel_rate=0.8):
	masks_index=[] #list of list of list, [mask[patch[cluster]]]
	used_channel=[]
	for i in range(len(dic_list)):
		print("Doing ", i)
		c_m0=dic_list[i]
		for c0 in list(c_m0.keys())[0:np.min([num_mask_each, len(c_m0)])]:
			if str(i)+"_"+str(c0) not in used_channel: #Start matching
				no_match=0
				masks_index_tmp=[]
				used_channel_tmp=[str(i)+"_"+str(c0)]
				for j in range(len(dic_list)):
					c_m1=dic_list[j]
					tmp1=[(k, np.max(np.abs(np.array(v)-np.array(c_m0[c0])))) for k, v in c_m1.items()]
					tmp2=list(filter(lambda x: x[1]<mapping_threshold1, tmp1))
					if len(tmp2)>0:
						masks_index_tmp.append([x[0] for x in tmp2])
						used_channel_tmp+=[str(j)+"_"+str(x[0])for x in tmp2]
						#print("Mapping: patch "+str(i)+" mask "+str(c0)+" to patch "+str(j)+" mask "+str(tmp[0][0])+"diff = "+str(np.round(tmp[0][1], 3)))
					else:
						tmp2=list(filter(lambda x: x[1]<mapping_threshold2, tmp1))
						if len(tmp2)>0:
							tmp2.sort(key=lambda x: x[1])
							masks_index_tmp.append([tmp2[0][0]])
							used_channel_tmp.append(str(j)+"_"+str(tmp2[0][0]))
						else:
							masks_index_tmp.append([])
							no_match+=1
				new_rate=1-len(set(used_channel_tmp)&set(used_channel))/len(used_channel_tmp)
				if no_match<(len(dic_list)*min_no_mach_rate) and (new_rate >= min_unused_channel_rate):
					print("Adding, ", "no match counts:",no_match, "new rate:",new_rate)
					masks_index.append(masks_index_tmp)
					used_channel+=used_channel_tmp
				else:
					print("Not adding, ", "no match counts:",no_match, "new rate:",new_rate)
	return masks_index

def Extract_Masks(masks_index, pred_file_locs, patch_size):
	num_masks=len(masks_index)
	masks=np.zeros([num_masks, len(masks_index[0]), patch_size, patch_size])
	for i in range(num_masks):
		print("Extracting mask ", i)
		mask_index=masks_index[i]
		for j in range(len(mask_index)):
			pred=np.load(pred_file_locs[j])
			mask = pred.reshape(patch_size, patch_size).astype(np.uint8)
			mask=1*np.isin(mask, mask_index[j])
			masks[i, j, :, :]=mask
	return masks

def Combine_Masks(masks, patch_info, img_size0, img_size1):
	#Combine masks to WSI
	patch_size=masks.shape[2]
	d0=int(np.ceil(img_size0/patch_size)*patch_size)
	d1=int(np.ceil(img_size1/patch_size)*patch_size)
	combined_masks=np.zeros([masks.shape[0], d0, d1])
	for i in range(masks.shape[0]): #Each mask
		print("Combining mask ", i)
		for j in range(masks.shape[1]): #Each patch
			info=patch_info.iloc[j]
			x, y=int(info["pixel_x"]), int(info["pixel_y"])
			combined_masks[i, int(x-patch_size/2):int(x+patch_size/2), int(y-patch_size/2):int(y+patch_size/2)]=masks[i, j, :, :]
	combined_masks=combined_masks[:, 0:img_size0, 0:img_size1]
	return combined_masks

def Extract_CC_Features_each_CC(labels):
	#Check all attr
	#[attr for attr in dir(region_props[0]) if not attr.startswith('__')]
	#labels size1 x size2
	pnames=["label","area",  "major_axis_length", "minor_axis_length",  "solidity"]
	region_props = regionprops(labels.astype(int))
	label_area=[(prop.label, prop.area) for prop in region_props]
	ret={}
	for name in pnames:
		ret[name]=[]
		for i in range(len(region_props)):
			ret[name].append(getattr(region_props[i],name))
	ret=pd.DataFrame(ret)
	return ret


# Convert to spot level¶
def extract_color(x_pixel=None, y_pixel=None, image=None, beta=49):
	beta_half=round(beta/2)
	g=[]
	for i in range(len(x_pixel)):
		max_x=image.shape[0]
		max_y=image.shape[1]
		nbs=image[max(0,x_pixel[i]-beta_half):min(max_x,x_pixel[i]+beta_half+1),max(0,y_pixel[i]-beta_half):min(max_y,y_pixel[i]+beta_half+1)]
		g.append(np.mean(nbs))
	c3=np.array(g)
	return c3

