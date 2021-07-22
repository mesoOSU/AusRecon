# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 09:25:10 2021

@author: arger
"""
import numpy as np


def dummy_header(header,Aus1_Mart2=True,sample_id='python_dummy_header',scan_id=00):
    if Aus1_Mart2 == True:
        P1 = 1
        P2 = 2
    else:
        P1=2
        P2=1
    grid     = header.split("GRID:")[1].split("#")[0].split('\n')[0]
    x_step   = header.split("XSTEP:")[1].split("#")[0].split('\n')[0]
    y_step   = header.split("YSTEP:")[1].split("#")[0].split('\n')[0]
    odd_col  = header.split("NCOLS_ODD:")[1].split("#")[0].split('\n')[0]
    even_col = header.split("NCOLS_EVEN:")[1].split("#")[0].split('\n')[0]
    rows     = header.split("NROWS:")[1].split("#")[0].split('\n')[0]
    header = """# TEM_PIXperUM          1.000000
# x-star                0.468457
# y-star                0.480084
# z-star                0.556985
# WorkingDistance       20.990000
#
# Phase {}
# MaterialName  	Austenite
# Formula     	Fe
# Info 		
# Symmetry              43
# LatticeConstants      3.650 3.650 3.650  90.000  90.000  90.000
# NumberFamilies        4
# hklFamilies   	 2  1  1 1 0.000000 1
# hklFamilies   	 3  1  0 1 0.000000 1
# hklFamilies   	 1  1  0 1 0.000000 1
# hklFamilies   	 2  0  0 1 0.000000 1
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# Categories0 0 0 0 0 
#
# Phase {}
# MaterialName  	Martensite
# Formula     	Fe
# Info 		
# Symmetry              43
# LatticeConstants      2.870 2.870 2.870  90.000  90.000  90.000
# hklFamilies   	 1  1  1 1 0.000000 1
# hklFamilies   	 2  0  0 1 0.000000 1
# hklFamilies   	 2  2  0 1 0.000000 1
# hklFamilies   	 3  1  1 1 0.000000 1
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# ElasticConstants 	-1.000000 -1.000000 -1.000000 -1.000000 -1.000000 -1.000000
# Categories0 0 0 0 0 
#
#
# GRID: {}
# XSTEP: {}
# YSTEP: {}
# NCOLS_ODD: {}
# NCOLS_EVEN: {}
# NROWS: {}
#
# OPERATOR: 	Austin Gerlt
#
# SAMPLEID: 	{}
#
# SCANID: 	{}
#""".format(P1,P2,grid,x_step,y_step,even_col,odd_col,rows,sample_id,scan_id)
    return(header)

        
with open('Map1.ang') as myfile:
    top = [next(myfile) for x in range(1000)]
    myfile.close()

original_header = ''.join(np.array(top)[np.array([x[0]=='#' for x in top])])
new_header = dummy_header(original_header)
dat = np.loadtxt('Map1.ang')
data  = dat[:,0:8] 
a = np.zeros([data.shape[0],6])
appendix = a+np.array([1, 0, 0, 0, 0, 0])
body = np.hstack([data,appendix])
fmt_str = '%0.5f %0.5f %0.5f %0.5f %0.5f %0.1f %0.1f %0.0f %0.0f %0.0f %0.0f %0.0f %0.0f %0.0f'

np.savetxt('test.ang',body,fmt=fmt_str,header=new_header,comments = "")