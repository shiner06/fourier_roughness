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
# NOTE: Axis of cylindrical body must be oriented along x axis

from initialize import mesh_reader, point_cloud_writer
from excel_writer import exporter
from bar_graph_3d import bar_graph
from example_plots import line_plot, contour_plot
import open3d as o3d
import os
import time
import numpy as np
import perturb

start_time = time.time()

input_path  = os.getcwd()+"/"
output_path = os.getcwd()+"/"
mesh_file   = "Common_Poly_Sharp_Tip_Rigid_10K_Rotated.STL"

##############################################################################
# Alter these fourier series function parameters
# The length of the sharp tip is 7.5 inches
# 7.5 = 2*N * lambda_k
N=20                # Upper limit of summation
M=4                 # Upper limit of summation
lambda_k=250        # Roughness design minimum wavelength in mils
##############################################################################

# Read an .STL or .PLY mesh, pass to point cloud writer
# mesh = mesh_reader(input_path, mesh_file)

#  Poisson disk sampling to create point cloud
# sample_points = 1000c

# pcd = point_cloud_writer(mesh, sample_points)

# Create amplitude and phase shift matrices
# to recompile perturb, enter: f2py -m perturb -c perturb.f90
print("Preparing randomized amplitudes and phases...")
A, phi = perturb.wave_prep(N, M, lambda_k)

# Plot a representative line plot
line_plot(phi, A, N, M, lambda_k)

# Plot a representative contour in cartesian coordinates
contour_plot(phi, A, N, M, lambda_k)

# Plot the A matrix
# bar_graph(A, N, M, lambda_k)

# Perturb point cloud with axisymmetric fourier series
# print("Perturbing point cloud...")
# pts = np.float32(np.asarray(pcd.points))
# pts_perturbed = perturb.axi_fourier_series(pts, phi, A, lambda_k, sample_points, N, M)

# Clear old elements
# pcd.points.clear()

# Replace with new perturbed elements
# pcd.points.extend(pts_perturbed)

# Write a new point cloud
# o3d.io.write_point_cloud("point_cloud_perturbed.xyz", pcd)

# print(pcd)

# Export fourier parameters to excel file
# exporter(A, phi, lambda_k, N, M, pts_perturbed, output_path)

# np.savetxt("coords.csv", pts_perturbed, delimiter=",")

print("Calculations completed in %.2f seconds" % (time.time() - start_time))
