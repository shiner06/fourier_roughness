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
from perturb import cart_fourier_series

def line_plot(phi, A, N, M, lambda_k):
    
    # Build a 1D plot for the simplest example

    title = '1D Example'

    sp = 1000 # sample_points
    x = np.linspace(-N*lambda_k, N*lambda_k,sp)
    y = []
    N_array = np.linspace(-N,N,2*N+1)
    A_array = []
    B_array = []
    # phi_array = np.zeros([2*N+1,])
    phi_array = lambda_k * np.random.randn((2*N+1),) / 4
    for n in N_array:
        A_array.append(lambda_k / 2 * n**2 / N**2)
    for xi in x:
        B_array = np.cos(2*np.pi*N_array*xi / (lambda_k * N) + phi_array)
        y.append((sum(A_array *  B_array)))
    
    fig,ax=plt.subplots(1,1)
    plt.plot(x,y,color='red')
    plt.grid(visible=True, which='Major', axis='both', linewidth=0.5)
    ax.set_xticks(     np.linspace(-N*lambda_k     , N*lambda_k     , 5))
    ax.set_xticklabels(np.linspace(-N*lambda_k/1000, N*lambda_k/1000, 5), fontname='Times New Roman', fontsize=24) 
    ax.set_yticks(     np.linspace(  -lambda_k     , lambda_k       , 5))
    ax.set_yticklabels(np.linspace(  -lambda_k     , lambda_k       , 5), fontname='Times New Roman', fontsize=24)
    ax.set_xlabel('x, inches'                                           , fontname='Times New Roman', fontsize=28)
    ax.set_ylabel('Roughness height, mils'                              , fontname='Times New Roman', fontsize=28)
    # ax.set_title(title                                                                         , fontname='Times New Roman', fontsize=20, fontweight='bold')
    
    # get current figure
    figure = plt.gcf()
    figure.set_size_inches(8, 6)
    
    # when saving, specify the DPI
    figure.savefig(title + '.png', bbox_inches="tight", dpi = 200)




def contour_plot(phi, A, N, M, lambda_k):

    # Contour Plot
    title = '2D Example'
    
    # Build a fake pts matrix
    samples = 512
    sample_points = samples**2
    xpos = np.linspace(-N*lambda_k, N*lambda_k, samples)
    ypos = np.linspace(-M*lambda_k, M*lambda_k, samples)
    xpos, ypos = np.meshgrid(xpos, ypos)
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    pts = np.zeros([sample_points, 2])
    pts[:,0] = xpos
    pts[:,1] = ypos
    
    x, y, z = cart_fourier_series(pts, phi, A, lambda_k, sample_points, N, M)
    
    xpos = np.reshape(x, (samples, samples))
    ypos = np.reshape(y, (samples, samples))
    zpos = np.reshape(z, (samples, samples))
    
    fig,ax=plt.subplots(1,1)
    cp =  ax.contourf(xpos, ypos, zpos, cmap=plt.cm.jet, levels=np.linspace(-lambda_k/2/1000-0.01, lambda_k/2/1000+0.01, 500))
    cb = fig.colorbar(cp)
    cb.set_ticks(      np.linspace(  -lambda_k/2/1000, lambda_k/2/1000, 3))
    cb.set_ticklabels( np.linspace(  -lambda_k/2     , lambda_k/2     , 3), fontname='Times New Roman', fontsize=24)
    ax.set_xticks(     np.linspace(-N*lambda_k/1000, N*lambda_k/1000, 5))
    ax.set_xticklabels(np.linspace(-N*lambda_k/1000, N*lambda_k/1000, 5), fontname='Times New Roman', fontsize=24)
    ax.set_yticks(     np.linspace(-M*lambda_k/1000, M*lambda_k/1000, 5))
    ax.set_yticklabels(np.linspace(-M*lambda_k/1000, M*lambda_k/1000, 5), fontname='Times New Roman', fontsize=24)
    cb.set_label('Roughness height, mils'                                         , fontname='Times New Roman', fontsize=28)
    ax.set_xlabel('x, inches'                                                     , fontname='Times New Roman', fontsize=28)
    ax.set_ylabel('y, inches'                                                     , fontname='Times New Roman', fontsize=28)
    # ax.set_title(title                                                                             , fontname='Times New Roman', fontsize=20, fontweight='bold')
    plt.gca().set_aspect('equal')
    
    # get current figure
    figure = plt.gcf() 
    figure.set_size_inches(8, 6)
    
    # when saving, specify the DPI
    figure.savefig(title + '.png', bbox_inches="tight", dpi = 200)
