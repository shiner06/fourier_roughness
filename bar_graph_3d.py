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

def bar_graph(A):
    fig = plt.figure()
    ax = Axes3D(fig)
    
    lx= len(A[0])            # Work out matrix dimensions
    ly= len(A[:,0])
    xpos = np.arange(0,lx,1)    # Set up a mesh of positions
    ypos = np.arange(0,ly,1)
    xpos, ypos = np.meshgrid(xpos+0.25, ypos+0.25)
    
    xpos = xpos.flatten()   # Convert positions to 1D array
    ypos = ypos.flatten()
    zpos = np.zeros(lx*ly)
    
    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = A.flatten()
    
    ax.bar3d(xpos,ypos,zpos, dx, dy, dz)
    
    plt.show()