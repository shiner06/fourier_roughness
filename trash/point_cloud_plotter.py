#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tues Aug 23 2022

@author: 
    
    Jack Shine
    
    Graduate Student Researcher            
    
    National Aerothermochemistry and Hypersonic Labratory       
    
    Texas A & M University  
"""
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math

def cloud_plot(pts, pts_perturbed):
        
    # Define nominal sampled points in cartesian
    x = pts[:,0]
    y = pts[:,1]
    z = pts[:,2]
    
    # Define perturbed points in cartesian
    xp = pts_perturbed[:,0]
    yp = pts_perturbed[:,1]
    zp = pts_perturbed[:,2]
    
    th = []
    r  = []
    # Map cartesian coordinates into cylindrical coordinate system
    for i in range(len(pts)):
      th.append(math.atan2(z[i], y[i]))
      r.append( np.sqrt(y[i]**2 + z[i]**2))
    
    th_perturbed = []
    r_perturbed  = []
    # Map cartesian coordinates into cylindrical coordinate system
    for i in range(len(pts_perturbed)):
      th_perturbed.append(math.atan2(zp[i], yp[i]))
      r_perturbed.append( np.sqrt(yp[i]**2 + zp[i]**2))
      
    # Creat a color assignment defining the extent of perturbation
    subtracted   = np.subtract(r_perturbed, r)
    absval       = np.abs(subtracted)
    fracs        = subtracted.astype(float)/absval.max()
    norm         = colors.Normalize(fracs.min(), fracs.max())
    color_values = cm.jet(norm(fracs.tolist()))
        
    # Contour Plot
    title = 'Perturbed Point Cloud'
    
    fig = plt.figure()
    ax  = Axes3D(fig)
    ax.set_xlabel('x (in)')
    ax.set_ylabel('y (in)')
    ax.set_zlabel('z (in)')
    ax.view_init(40, 245)
    ax.set_title(title)
    ax.scatter(x, y, z, s=20, depthshade=True, c = color_values)
    ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
    
    # get current figure
    figure = plt.gcf() 
    figure.set_size_inches(8, 6)
    
    # when saving, specify the DPI
    figure.savefig(title + '.png', dpi = 200)
    