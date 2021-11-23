# -*- coding: utf-8 -*-
"""
Created on Mon Sep 20 16:34:20 2021

pythonic graph cut algorithm (single layer, node/source example)

@author: arger
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
import networkx as nx
#from networkx import DiGraph
np.seterr(divide='ignore')

#clear out any open plots
plt.close('all')

# read in image that requires segmentation
#img = plt.imread("slice001.tiff")
img =np.mean( plt.imread("LabHEDM_zoom1.tif"),axis = 2)
#img = plt.imread("office.jpg")[500::5,100::5,2]
if img.max()>1:
    img = img/256
img_vec = img.reshape([-1,1])

# fit a 2-part Gaussian Mixture function to it (black matrix and white precipitate)
GMD = GaussianMixture(2).fit(img.reshape(-1,1))

#define the 2 gaussien functions from the gmd data
precip_fit = norm(loc = GMD.means_[1],scale = GMD.covariances_[1,0,0]**0.5)
matrix_fit = norm(loc = GMD.means_[0],scale = GMD.covariances_[0,0,0]**0.5)

#Plot the curve fits, just to see them
x_bar = np.arange(0,0.6,1/512)
x_bar = np.arange(0,1.0,1/512)
precip_curve =  precip_fit.pdf(x_bar)*GMD.weights_[1]
matrix_curve =  matrix_fit.pdf(x_bar)*GMD.weights_[0]
histo = plt.hist(img_vec, bins = x_bar,density = 'norm')
plt.plot(x_bar,precip_curve*2)
plt.plot(x_bar,matrix_curve*2)

#Assign each pixel a probablity that it belongs to either the matrix or source


# Now find the inverse deltas between colors for each in-plane connection
ID = np.arange(img.size).reshape(img.shape)
# next two objects create a 3-tuple of A,B, and 1/delta(AB). This gives similar
# neighbors a strong bond and dissimilar ones a weak bond. The third line 
# stacks them.
LR_edges = ID[:,:-1], ID[:,1:], (1/np.abs(img[:,:-1]-img[:,1:])).astype(int)
UD_edges = ID[:-1,:], ID[1:,:], (1/np.abs(img[:-1,:]-img[1:,:])).astype(int)
IP_edge_weights = np.vstack([np.vstack([x.flatten() for x in LR_edges]).T,
                             np.vstack([x.flatten() for x in UD_edges]).T])
# the code above gives identical neighbors a strong negative bond. (intentional
# artifact of casting "inf" as an int) this line replaces those values 500
IP_edge_weights[(IP_edge_weights[:,2]<0),2] = 500

# =======================
#Start making the graph.
DG = nx.Graph()
#add N+2 nodes for In-plane(0 to n-1), source(n) and sink(n+1) nodes  
DG.add_node(img.size+2)

#Alright, there is a DEFINITELY a more pythonic/efficient way to do this, but
#adding edges en-masse is weird, so here is a slowish "for" loop
IP_weight = 0
IP_scale = 1 
for entry in IP_edge_weights:
    DG.add_edge(entry[0],entry[1],capacity=float((entry[2]*IP_scale)+IP_weight))

# Make pixel-to-source weights from pdf
particle_weight =0
particle_scale = 25
particle_2_pix = ((precip_fit.pdf(img_vec)*particle_scale)+particle_weight).flatten()
    
# Make pixel-to-sink weights from pdf
matrix_weight = 0
matrix_scale = 250
matrix_2_pix = ((matrix_fit.pdf(img_vec)*matrix_scale)+matrix_weight).flatten()


# add those weights in
for i in np.arange(img_vec.size):
    DG.add_edge(img_vec.size+1,i,capacity=particle_2_pix[i])
    DG.add_edge(img_vec.size+2,i,capacity=matrix_2_pix[i])

cut, part = nx.minimum_cut(DG,img_vec.size+1,img_vec.size+2)

#particle_2_pix = ((1-img_vec)*particle_scale)+particle_weight
#matrix_2_pix = ((img_vec)*matrix_scale)+matrix_weight

# add those weights in
for i in np.arange(img_vec.size):
    DG.add_edge(img_vec.size+1,i,capacity=particle_2_pix[i])
    DG.add_edge(img_vec.size+2,i,capacity=matrix_2_pix[i])

cut, part = nx.minimum_cut(DG,img_vec.size+1,img_vec.size+2)

src = img_vec*1
snk = img_vec*1

snk_IDs = np.array(list(part[0]))
src_IDs = np.array(list(part[1]))
snk_IDs = snk_IDs[snk_IDs<=img_vec.size]
src_IDs = src_IDs[src_IDs<=img_vec.size]
src[src_IDs] = 0
snk[snk_IDs] = 0


plt.figure()
plt.imshow(img_vec.reshape(img.shape))#,cmap = 'gray')
plt.figure()
plt.imshow(src.reshape(img.shape))#,cmap = 'gray',vmax = 0.01)
plt.figure()
plt.imshow(snk.reshape(img.shape))#,cmap = 'gray',vmax = 0.01)


plt.figure()
plt.imshow(src.reshape(img.shape))#,cmap = 'gray')
plt.figure()
plt.imshow(snk.reshape(img.shape))#,cmap = 'gray')

