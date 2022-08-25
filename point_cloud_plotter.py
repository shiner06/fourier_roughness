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
import matplotlib.pyplot as plt
import numpy as np
import perturb

def point_cloud_plotter(pts, pts_perturbed):

    # Contour Plot
    title = 'Perturbed Point Cloud'
    
    # Map cartesian coordinates into cylindrical coordinate system
    # for i in sample_points-1
    #   th(i) = atan2(z(i),y(i))
    #   r(i) = sqrt(y(i)**2 + z(i)**2)
    # end do

    # # Build a fake pts matrix
    # samples = 512
    # sample_points = samples**2
    # xpos = np.linspace(-N*lambda_k, N*lambda_k, samples)
    # ypos = np.linspace(-M*lambda_k, M*lambda_k, samples)
    # xpos, ypos = np.meshgrid(xpos, ypos)
    # xpos = xpos.flatten()
    # ypos = ypos.flatten()
    # pts = np.zeros([sample_points, 2])
    # pts[:,0] = xpos
    # pts[:,1] = ypos
    
    # x, y, z = perturb.cart_fourier_series(pts, phi, A, lambda_k, sample_points, N, M)
    
    # xpos = np.reshape(x, (samples, samples))
    # ypos = np.reshape(y, (samples, samples))
    # zpos = np.reshape(z, (samples, samples))
    
    # fig,ax=plt.subplots(1,1)
    # cp =  ax.contourf(xpos, ypos, zpos, cmap=plt.cm.jet, levels=np.linspace(-lambda_k*1.25/1000, lambda_k*1.25/1000, 500))
    # cb = fig.colorbar(cp)
    # cb.set_ticks(                np.linspace(  -lambda_k/1000,   lambda_k/1000, 5),              fontname='Times New Roman', fontsize=16)
    # cb.set_ticklabels( np.around(np.linspace(  -lambda_k/1000,   lambda_k/1000, 5), decimals=3), fontname='Times New Roman', fontsize=16, fontweight='bold')
    # ax.set_xticks(               np.linspace(-N*lambda_k/1000, N*lambda_k/1000, 5)             , fontname='Times New Roman', fontsize=20)
    # ax.set_yticks(               np.linspace(-M*lambda_k/1000, M*lambda_k/1000, 5)             , fontname='Times New Roman', fontsize=16)
    # ax.set_xticklabels(np.around(np.linspace(-N*lambda_k/1000, N*lambda_k/1000, 5), decimals=1), fontname='Times New Roman', fontsize=16)
    # ax.set_yticklabels(np.around(np.linspace(-M*lambda_k/1000, M*lambda_k/1000, 5), decimals=1), fontname='Times New Roman', fontsize=16)
    # cb.set_label('Roughness height [in]'                                                       , fontname='Times New Roman', fontsize=18)
    # ax.set_xlabel('x (in)'                                                                     , fontname='Times New Roman', fontsize=18)
    # ax.set_ylabel('y (in)'                                                                     , fontname='Times New Roman', fontsize=18)
    # ax.set_title(title                                                                         , fontname='Times New Roman', fontsize=20, fontweight='bold')
    # plt.gca().set_aspect('equal')
    
    # # get current figure
    # figure = plt.gcf() 
    # figure.set_size_inches(8, 6)
    
    # # when saving, specify the DPI
    # figure.savefig(title + '.png', dpi = 200)