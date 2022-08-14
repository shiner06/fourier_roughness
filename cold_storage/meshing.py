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

def poisson_mesh(pcd, output_file):
    """
    Using a point cloud object (pcd), reconstruct a surface with the
    poisson surface reconstruction algorithm
    """
    
    print("Creating Surface Mesh...")
    
    # Create a new mesh
    mesh, densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=11)
    
    print("Removing spurious clusters...")
    with o3d.utility.VerbosityContextManager(o3d.utility.VerbosityLevel.Debug):
        triangle_clusters, cluster_n_triangles, cluster_area = (mesh.cluster_connected_triangles())
        triangle_clusters = np.asarray(triangle_clusters)
        cluster_n_triangles = np.asarray(cluster_n_triangles)
        largest_cluster_idx = cluster_n_triangles.argmax()
        triangles_to_remove = triangle_clusters != largest_cluster_idx
        mesh.remove_triangles_by_mask(triangles_to_remove)
        print("Largest cluster identified")
    
    # Remove any non-manifold edges
    if not mesh.is_edge_manifold():
        mesh.remove_non_manifold_edges()
        print("Non-manifold edges removed")
    
    # Subdivide with loop algorithm
    print("Subdividing for high resolution...")
    mesh = mesh.subdivide_loop(number_of_iterations=3)
    
    # Remove more non-manifold edges
    if not mesh.is_edge_manifold():
        mesh.remove_non_manifold_edges()
        print("More Non-manifold edges removed")
    
    # Calculate normals
    mesh.compute_vertex_normals()
    
    # Colorize the mesh
    mesh.paint_uniform_color([1, 1, 1])
    
    # Write a mesh file
    o3d.io.write_triangle_mesh(output_file, mesh)
    print("Surface mesh calculations complete")
    print(mesh)
    return mesh