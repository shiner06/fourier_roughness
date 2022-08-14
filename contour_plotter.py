#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 2022

@author: 
    
    Jack Shine
    
    Graduate Student Researcher            
    
    National Aerothermochemistry and Hypersonic Labratory       
    
    Texas A & M University  
"""
import matplotlib.pyplot as plt
import numpy as np
import perturb

def contour(phi, A, N, M, lambda_k):
    
    # X/Y Domain should run from -lambda_k*N to lambda_k * N
    # In other words, there are a maximum of 2*N modes each 
    # of length lambda_k running in one domain direction

    # Build a fake pts matrix
    samples = 512
    sample_grid = samples**2
    xpos = np.linspace(-N*lambda_k, N*lambda_k, samples)
    ypos = np.linspace(-M*lambda_k, M*lambda_k, samples)
    xpos, ypos = np.meshgrid(xpos, ypos)
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    pts = np.zeros([sample_grid, 2])
    pts[:,0] = xpos
    pts[:,1] = ypos

    x, y, z = perturb.cart_fourier_series(pts, phi, A, lambda_k, sample_grid, N, M)
    
    xpos = np.reshape(x, (samples, samples))
    ypos = np.reshape(y, (samples, samples))
    zpos = np.reshape(z, (samples, samples))
    
    fig,ax=plt.subplots(1,1)
    cp = ax.contourf(xpos, ypos, zpos, cmap=plt.cm.jet, levels=100)
    cb = fig.colorbar(cp)
    cb.set_label('Roughness height [in]')
    ax.set_title('Representative Roughness Contour Plot')
    ax.set_xlabel('x (in)')
    ax.set_ylabel('y (in)')
    plt.gca().set_aspect('equal')
    
    # get current figure
    figure = plt.gcf() 
    figure.set_size_inches(8, 6)
    
    # when saving, specify the DPI
    figure.savefig("representative_contour.png", dpi = 200)
