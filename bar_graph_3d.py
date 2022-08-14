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
    lx          = len(A[0])            
    ly          = len(A[:,0])
    xpos        = np.arange(0,lx,1)
    ypos        = np.arange(0,ly,1)
    xpos, ypos  = np.meshgrid(xpos, ypos)
    
    # Convert positions to 1D array
    xpos = xpos.flatten()   
    ypos = ypos.flatten()
    zpos = np.zeros(lx*ly)
    
    # Control the thickness/spacing of bars and tick marks
    scale_factor    = 0.5
    dx   = scale_factor * np.ones_like(zpos)
    dy   = dx.copy()
    dz   = A.flatten()
    
    offset          = dz + np.abs(dz.min())
    fracs           = offset.astype(float)/offset.max()
    norm            = colors.Normalize(fracs.min(), fracs.max())
    color_values    = cm.jet(norm(fracs.tolist()))

    ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=color_values)
    
    x_ticks = np.linspace(1,len(range(2*N+1)),int(2*N*scale_factor+1))
    y_ticks = np.linspace(1,len(range(2*M+1)),int(2*M*scale_factor+1))
    z_ticks = np.linspace(0,int(lambda_k),int(lambda_k/50)+1).astype(int)
    x_tick_labels = np.round(np.linspace(-N,N,int((2*N)*scale_factor)+1)).astype(int)
    y_tick_labels = np.round(np.linspace(-M,M,int((2*M)*scale_factor)+1)).astype(int)
    
    ax.w_xaxis.set_ticks(ticks = x_ticks, labels=x_tick_labels)
    ax.w_yaxis.set_ticks(ticks = y_ticks, labels=y_tick_labels)
    ax.set_zticks(z_ticks)
    ax.set_xlabel('N')
    ax.set_ylabel('M')
    ax.set_zlabel('Amplitude Coefficient (mils)')
    ax.view_init(50, -45)
    # plt.show()
    
    # get current figure
    figure = plt.gcf() 
    figure.set_size_inches(8, 6)
    
    # when saving, specify the DPI
    figure.savefig("amplitude_coefficients.png", dpi = 200)
