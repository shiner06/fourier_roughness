#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 15:00:10 2020

@author: 
    
    Jack Shine
    
    Graduate Student Researcher            
    
    National Aerothermochemistry and Hypersonic Labratory       
    
    Texas A & M University  
"""
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm

def bar_graph(A, N, M, lambda_k):
    fig = plt.figure()
    ax  = Axes3D(fig)
    
    # Work out matrix dimensions and create mesh
    lx          = np.size(A[0,:])
    ly          = np.size(A[:,0])
    xpos        = np.arange(0,lx,1)
    ypos        = np.arange(0,ly,1)
    xpos, ypos  = np.meshgrid(xpos, ypos)
    
    # Convert positions to 1D array
    xpos = xpos.flatten()   
    ypos = ypos.flatten()
    zpos = np.zeros(lx*ly)
    
    # Control the thickness/spacing of bars and tick marks
    scale_factor = 0.5
    dx           = scale_factor * np.ones_like(zpos)
    dy           = dx.copy()
    dz           = A.flatten()
    
    fracs           = dz.astype(np.float)/dz.max()
    norm            = colors.Normalize(fracs.min(), fracs.max())
    color_values    = cm.jet(norm(fracs.tolist()))

    ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=color_values)
    
    x_ticks = np.linspace(1,2*N+1,np.int(N/2 + 1))
    y_ticks = np.linspace(1,2*M+1,np.int(M/2 + 1))
    z_ticks = np.linspace(0,np.int(lambda_k),np.int(lambda_k/50)+1).astype(np.int)
    x_tick_labels = np.round(np.linspace(-N,N,np.int(N/2 + 1))).astype(np.int)
    y_tick_labels = np.round(np.linspace(-M,M,np.int(N/2 + 1))).astype(np.int)
    
    ax.w_xaxis.set_ticks(ticks = x_ticks, labels=x_tick_labels)
    ax.w_yaxis.set_ticks(ticks = y_ticks, labels=y_tick_labels)
    ax.set_zticks(z_ticks)
    ax.set_xlabel('N')
    ax.set_ylabel('M')
    ax.set_zlabel('Amplitude Coefficient (mils)')
    ax.view_init(45, -45)
    ax.dist = 12
    
    # get current figure
    fig.set_size_inches(8, 6)
    
    # when saving, specify the DPI
    fig.savefig("Amplitude Coefficients.png", dpi = 200)
