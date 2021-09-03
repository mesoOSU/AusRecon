# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 02:42:15 2021

@author: arger
"""
import numpy as np
from scipy.spatial.transform import Rotation as R
from scipy.io import loadmat

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection

def R_to_IPF(Rot,miller = np.array([0,0,1])):
    # two quick notes: 
        #1) Only works for cubic (Other systems require reciprical lattice)
        #2) using cheat eric taught me where i reorder xyz and magic happens
    true_xyz = np.abs(np.matmul(miller,Rot.as_matrix()))
    xyz = np.sort(true_xyz,axis = 1)
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]
    x_prime = x/(np.abs(z)+1)
    y_prime = y/(np.abs(z)+1)
    lo = np.min(np.vstack([x_prime,y_prime]),axis = 0)
    hi = np.max(np.vstack([x_prime,y_prime]),axis = 0)
    reg = np.abs(z)
    rgb = np.vstack([z-y,y-x,x]).T
    rgb = rgb/np.max(rgb,axis=0)
    rgb = (rgb.T/np.max(rgb,axis=1)).T
    return(rgb)

def xyz_to_xy_prime(x,y,z):
    x_prime = xyz[0]/(1+xyz[2])
    y_prime = xyz[1]/(1+xyz[2])
    return(x_prime,y_prime)

def xyz_to_rgb(x,y,z):
    x_prime = xyz[0]/(1+xyz[2])
    y_prime = xyz[1]/(1+xyz[2])
    return(x_prime,y_prime)

# mat = loadmat("EBSD/AF_001/RAF_001_HW10_IP20-40_OP10-100.mat")


plt.close ('all')

fig,ax = plt.subplots()
thetas = np.linspace(0,2*np.pi,1000)
x = np.cos(thetas)
y = np.sin(thetas)
z = y*0

ax.plot(x,y,c = 'k')
#ax.plot(x/2,y/2,c = 'k')

Rx = R.from_euler('zxz',[0 ,45,0],degrees = True).as_matrix()
Ry = R.from_euler('zxz',[0,45,90],degrees = True).as_matrix()
xyz = np.matmul(Rx,np.vstack([x,y,z]))
x_prime = xyz[0]/(np.abs(xyz[2])+1)
y_prime = xyz[1]/(np.abs(xyz[2])+1)
ax.plot(x_prime,y_prime,c = 'k')
xyz = np.matmul(Ry,np.vstack([x,y,z]))
x_prime = xyz[0]/(np.abs(xyz[2])+1)
y_prime = xyz[1]/(np.abs(xyz[2])+1)
ax.plot(x_prime,y_prime,c = 'k')


sr12 = 1/np.sqrt(2)
ax.plot([0,0],[-1,1],c = 'k')
ax.plot([-1,1],[0,0],c = 'k')
ax.plot([-sr12,sr12],[-sr12,sr12],c = 'k')
ax.plot([sr12,-sr12],[-sr12,sr12],c = 'k')

rotations = R.from_quat(1-np.random.rand(100000,4)*2)
xyz = np.matmul(np.array([0,0,1]),rotations.as_matrix()).T
x_prime = xyz[0]/(np.abs(xyz[2])+1)
y_prime = xyz[1]/(np.abs(xyz[2])+1)

#abc = np.sort(np.abs(np.vstack([x_prime,y_prime,xyz[2]])),axis =0)#/np.sum(abc, axis = 0)
abc = np.sort(np.abs(xyz),axis =0)#/np.sum(abc, axis = 0)
rgb = np.array([abc[2]-abc[1],abc[1]-abc[0],abc[0]]).T
rgb = (rgb.T/np.max(rgb,axis=1)).T

plt.scatter(x_prime,y_prime, c = rgb,s = 1)
ax.set_aspect(1)



t = (xyz[0]>0)*(xyz[1]>0)*(xyz[2]>0)*(xyz[0]>xyz[1])*(xyz[1]<(xyz[2]+1))*(xyz[0]<(xyz[2]+1))

plt.figure()
plt.scatter(x_prime[t],y_prime[t], c = rgb[t],s = 1)
ax.set_aspect(1)
