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



def poisson_mesh(pcd, pcd_perturbed):
    """
    Using the perturbed point cloud object (pcd_perturbed), reconstruct a surface with the
    poisson surface reconstruction algorithm
    """
    
    ###########################################################################
    # Mesh the Perturbed Point Cloud
    ###########################################################################
    
    print("Creating perturbed surface mesh...")
    
    # Create a new mesh for the perturbed point cloud
    perturbed_mesh, perturbed_mesh_densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd_perturbed, depth=10)
    
    # with o3d.utility.VerbosityContextManager(o3d.utility.VerbosityLevel.Debug):
    triangle_clusters, cluster_n_triangles, cluster_area = (perturbed_mesh.cluster_connected_triangles())
    triangle_clusters = np.asarray(triangle_clusters)
    cluster_n_triangles = np.asarray(cluster_n_triangles)
    largest_cluster_idx = cluster_n_triangles.argmax()
    if triangle_clusters.any() != largest_cluster_idx:
        print("Removing spurious clusters...")
        triangles_to_remove = triangle_clusters != largest_cluster_idx
        perturbed_mesh.remove_triangles_by_mask(triangles_to_remove)
        print("Largest cluster maintained")
    
    # Remove non-Manifold edges if present
    if not perturbed_mesh.is_edge_manifold():
        perturbed_mesh.remove_non_manifold_edges()
        print("Non-manifold edges removed from perturbed mesh")
    
    # Subdivide with loop algorithm
    print("Subdividing perturbed mesh for high resolution...")
    perturbed_mesh = perturbed_mesh.subdivide_loop(number_of_iterations=2)
    
    # Remove non-manifold edges if present after subdividing
    if not perturbed_mesh.is_edge_manifold():
        perturbed_mesh.remove_non_manifold_edges()
        print("Non-manifold edges removed from perturbed mesh after subdivision")
    
    # Calculate normals
    perturbed_mesh.compute_vertex_normals()
    
    pts_perturbed = np.asarray(perturbed_mesh.vertices)

    print("There are " + str(len(pts_perturbed)) + " vertices in the perturbed mesh")             



    ###########################################################################
    # Mesh the Nominal Point Cloud
    ###########################################################################
        
    # Create a new mesh for the perturbed point cloud
    mesh, mesh_densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=10)
    
    # Remove non-Manifold edges if present
    if not mesh.is_edge_manifold():
        mesh.remove_non_manifold_edges()
        print("Non-manifold edges removed from smooth mesh")
        
    # Subdivide with loop algorithm
    print("Subdividing smooth mesh for high resolution...")
    mesh = mesh.subdivide_loop(number_of_iterations=1)
    
    # Remove non-manifold edges if present after subdividing
    if not mesh.is_edge_manifold():
        mesh.remove_non_manifold_edges()
        print("Non-manifold edges removed from smooth mesh after subdivision")
    
    # Calculate normals
    mesh.compute_vertex_normals()
    mesh.paint_uniform_color([1,1,1])
    
    pts = np.asarray(mesh.vertices)
    
    print("There are " + str(len(pts)) + " vertices in the smooth mesh") 
    
    
    
    ###########################################################################
    # Colorize the mesh    
    ###########################################################################
    
    # Find the nearest neighbor for each point in pts_perturbed:
    kdtree = KDTree(pts)
    dist, idx = kdtree.query(pts_perturbed,1)
    
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
    for i, pt in enumerate(pts):
      th.append(math.atan2(z[i], y[i]))
      r.append( np.sqrt(y[i]**2 + z[i]**2))
    
    th_perturbed = []
    r_perturbed  = []
    # Map cartesian coordinates into cylindrical coordinate system
    for i, pt_perturbed in enumerate(pts_perturbed):
      th_perturbed.append(math.atan2(zp[i], yp[i]))
      r_perturbed.append( np.sqrt(yp[i]**2 + zp[i]**2))
    
    # Uncomment this section if no longer displacing the pattern by lambda_k/2
    # If the radius of the perturbed point is less than the nominal radius of 
    # its nearest neighbor, multiply the distance by -1, otherwise do nothing
    # for i, each_perturbed_radius in enumerate(r_perturbed):
    #     if each_perturbed_radius < r[idx[i]]:
    #         dist[i] = dist[i] * -1
    
    # Create a color assignment defining the extent of perturbation
    # print(np.amax(np.abs(dist)))
    fracs        = dist / np.amax(np.abs(dist))
    norm         = colors.Normalize(fracs.min(), fracs.max())
    mesh_colors  = color_map.jet(norm(fracs.tolist()))
    mesh_colors  = mesh_colors[:, :3]
    perturbed_mesh.vertex_colors = o3d.utility.Vector3dVector(mesh_colors)
    
    
    
    ###########################################################################
    # Save Images and Mesh STL
    ###########################################################################
    
    # Save the Perturbed Mesh Contour Image
    vis = o3d.visualization.VisualizerWithEditing()
    vis.create_window(visible=False)
    vis.add_geometry(perturbed_mesh)
    vis.update_geometry(perturbed_mesh)
    view = vis.get_view_control()
    render = vis.get_render_option()
    render.background_color = [0.4, 0.4, 0.4]
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
    vis.add_geometry(mesh)
    vis.update_geometry(mesh)
    view = vis.get_view_control()
    render = vis.get_render_option()
    render.background_color = [0.4, 0.4, 0.4]
    vis.reset_view_point(True)
    view.set_zoom(0.5)
    vis.update_renderer()
    vis.poll_events()
    vis.capture_screen_image(os.getcwd() + '/Smooth Mesh.png')
    vis.destroy_window()
        
    # Write Mesh Files
    o3d.io.write_triangle_mesh('Perturbed Mesh.stl', perturbed_mesh)
    o3d.io.write_triangle_mesh('Smooth Mesh.stl'  , mesh)
    print("Perturbed surface meshing complete")
    print(perturbed_mesh)