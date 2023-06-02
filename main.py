#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  8 15:21:43 2020

@author:

    Jack Shine

    Graduate Student Researcher

    National Aerothermochemistry and Hypersonic Labratory

    Texas A & M University
"""
from perturb import wave_prep, axi_fourier_series
from initialize import mesh_reader, point_cloud_writer
from excel_writer import exporter
from bar_graph_3d import bar_graph
from example_plots import line_plot, contour_plot
# from meshing import poisson_mesh
import open3d as o3d
import os
import time
import numpy as np
import copy

start_time = time.time()

input_path  = os.getcwd()+"/"
output_path = os.getcwd()+"/"
mesh_file   = "OgiveForebody_OI.STL"
perturbed_mesh_file = "OgiveForebody_OI_Perturbed.STL"

##############################################################################
# Alter these fourier series function parameters
# The length of the sharp tip is 7.5", the base circumference is 6.855", roughly 7" inch average
# 3.5 / 0.035 = 100, 3.5 / 0.050 = 70, and 3.5 / 0.070 = 50
N=4              # Upper limit of summation
M=4              # Upper limit of summation
lambda_k=250     # The desired roughness wavelength (trough to trough) diameter in mils
##############################################################################

# Read an .STL or .PLY mesh, pass to point cloud writer
# mesh = mesh_reader(input_path, mesh_file)

#  Poisson disk sampling to create point cloud
# sample_points = 10000000

# pcd = point_cloud_writer(mesh, sample_points)

# Create amplitude and phase shift matrices
# to recompile perturb, enter: f2py -m perturb -c perturb.f90
A, phi = wave_prep(N, M, lambda_k)

# Plot a representative line plot
line_plot(phi, A, N, M, lambda_k)

# Plot a representative contour in cartesian coordinates
contour_plot(phi, A, N, M, lambda_k)

# Plot the A matrix
# bar_graph(A, N, M, lambda_k)

# Perturb point cloud with axisymmetric fourier series
# pts = np.asarray(pcd.points)
# pts_perturbed = axi_fourier_series(pts, phi, A, lambda_k, sample_points, N, M)

# Shallow copy pcd to pcd_perturbed to create new object without affecting pts,
# then clear elements
# pcd_perturbed = copy.copy(pcd)
# pcd_perturbed.points.clear()

# Replace with new perturbed elements
# pcd_perturbed.points.extend(pts_perturbed)

# Write a new point cloud
# o3d.io.write_point_cloud("point_cloud_perturbed.xyz", pcd_perturbed)

# Export fourier parameters to excel file
# exporter(A, phi, lambda_k, N, M, output_path)

# Export the perturbed coordinates
# np.savetxt("coords.csv", pts_perturbed, delimiter=",")

print("Calculations completed in %.2f seconds" % (time.time() - start_time))
