#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 12:30:30 2020

@author: 
    
    Jack Shine
    
    Graduate Student Researcher            
    
    National Aerothermochemistry and Hypersonic Labratory       
    
    Texas A & M University  
"""
import open3d as o3d
import numpy as np
import matplotlib.cm as color_map
import math
import matplotlib.colors as colors
import os
from scipy.spatial import KDTree

input_path          = os.getcwd()+"/"
output_path         = os.getcwd()+"/"
smooth_mesh_file    = "Common Poly Medium Tip.STL"
perturbed_mesh_file = "2.0x Perturbed Mesh.STL"

###########################################################################
# Import and color the smooth mesh
###########################################################################

# Create a new mesh for the perturbed point cloud
smooth_mesh = o3d.io.read_triangle_mesh(input_path+smooth_mesh_file)

# Sample a new point cloud from mesh
print("Sampling smooth mesh to create point cloud")
pcd = smooth_mesh.sample_points_uniformly(number_of_points=10000000)

# Estimate normals
# pcd.estimate_normals()
# pcd.orient_normals_consistent_tangent_plane(100)

# Create a new mesh for the perturbed point cloud
print("Creating new smooth mesh from sampled point cloud")
smooth_mesh, mesh_densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=10)

# Remove non-Manifold edges if present
if not smooth_mesh.is_edge_manifold():
    smooth_mesh.remove_non_manifold_edges()
    print("Non-manifold edges removed from smooth mesh")
    
# Subdivide with loop algorithm
print("Subdividing smooth mesh for high resolution...")
smooth_mesh = smooth_mesh.subdivide_loop(number_of_iterations=2)

# Remove non-manifold edges if present after subdividing
if not smooth_mesh.is_edge_manifold():
    smooth_mesh.remove_non_manifold_edges()
    print("Non-manifold edges removed from smooth mesh after subdivision")

# Compute normals for better rendering
print("Estimating smooth mesh normals")
smooth_mesh.compute_vertex_normals()

# Paint the smooth mesh
smooth_mesh.paint_uniform_color([0,0,1.0])

# Get vertex coordinates
pts = np.asarray(smooth_mesh.vertices)

print("There are " + str(len(pts)) + " vertices in the smooth mesh") 



###########################################################################
# Import and color the perturbed mesh
###########################################################################

# Create a new mesh for the perturbed point cloud
perturbed_mesh = o3d.io.read_triangle_mesh(input_path+perturbed_mesh_file)

# Compute normals for better rendering
print("Estimating perturbed mesh normals")
perturbed_mesh.compute_vertex_normals()

# Get vertex coordinates
pts_perturbed = np.asarray(perturbed_mesh.vertices)

print("There are " + str(len(pts_perturbed)) + " vertices in the perturbed mesh") 

# Find the nearest neighbor for each point in pts_perturbed:
kdtree = KDTree(pts)
dist, idx = kdtree.query(pts_perturbed,1)

# Create a color assignment defining the extent of perturbation
# print(np.amax(np.abs(dist)))
fracs        = dist / np.amax(np.abs(dist))
norm         = colors.Normalize(fracs.min(), fracs.max())
vert_colors  = color_map.jet(norm(fracs.tolist()))
vert_colors  = vert_colors[:, :3]
perturbed_mesh.vertex_colors = o3d.utility.Vector3dVector(vert_colors)



###########################################################################
# Save Images
###########################################################################

# Save the Perturbed Mesh Contour Image
vis = o3d.visualization.VisualizerWithEditing()
vis.create_window(visible=False)
vis.add_geometry(perturbed_mesh)
vis.update_geometry(perturbed_mesh)
view = vis.get_view_control()
render = vis.get_render_option()
render.background_color = [1.0, 1.0, 1.0]
vis.reset_view_point(True)
view.set_zoom(0.5)
vis.update_renderer()
vis.poll_events()
vis.capture_screen_image(os.getcwd() + '/Perturbed Mesh.png')
vis.destroy_window()
vis = []
view = []
render = []

# Save the Smooth Mesh Contour Image
vis = o3d.visualization.VisualizerWithEditing()
vis.create_window(visible=False)
vis.add_geometry(smooth_mesh)
vis.update_geometry(smooth_mesh)
view = vis.get_view_control()
render = vis.get_render_option()
render.background_color = [1.0, 1.0, 1.0]
vis.reset_view_point(True)
view.set_zoom(0.5)
vis.update_renderer()
vis.poll_events()
vis.capture_screen_image(os.getcwd() + '/Smooth Mesh.png')
vis.destroy_window()

# o3d.io.write_triangle_mesh('Smooth Mesh.STL'  , mesh)
print("Perturbed mesh visualization complete")
print(perturbed_mesh)