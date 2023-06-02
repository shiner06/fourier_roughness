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
    scale_factor = 0.4
    dx           = scale_factor * np.ones_like(zpos)
    dy           = dx.copy()
    dz           = A.flatten()
    
    fracs        = dz.astype(np.float)/dz.max()
    norm         = colors.Normalize(fracs.min(), fracs.max())
    color_values = cm.jet(norm(fracs.tolist()))

    ax.bar3d(xpos,ypos,zpos, dx, dy, dz, color=color_values, alpha=1)
    
    x_ticks = np.linspace(1,2*N+1,np.int(9))
    y_ticks = np.linspace(1,2*M+1,np.int(9))
    z_ticks = np.linspace(0,np.int(lambda_k/2),np.int(lambda_k/50)+1).astype(np.int)
    x_tick_labels = np.round(np.linspace(-N,N,np.int(9))).astype(np.int)
    y_tick_labels = np.round(np.linspace(-M,M,np.int(9))).astype(np.int)
    
    ax.w_xaxis.set_ticks(ticks = x_ticks, labels = x_tick_labels, fontname='Times New Roman', fontsize=20)
    ax.w_yaxis.set_ticks(ticks = y_ticks, labels = y_tick_labels, fontname='Times New Roman', fontsize=20)
    ax.w_zaxis.set_ticks(ticks = z_ticks, labels = z_ticks, fontname='Times New Roman', fontsize=20)
    ax.set_xlabel('N', fontname='Times New Roman', fontsize=20, labelpad = 20)
    ax.set_ylabel('M', fontname='Times New Roman', fontsize=20, labelpad = 20)
    ax.set_zlabel('\n\nAmplitude \nCoefficients, mils', fontname='Times New Roman', fontsize=20)
    ax.view_init(55, -44)
    ax.dist = 12
    
    # get current figure
    fig.set_size_inches(8, 6)
    
    # when saving, specify the DPI
    fig.savefig("Amplitude Coefficients.png", dpi = 200)
