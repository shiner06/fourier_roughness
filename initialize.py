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

import open3d as o3d

def mesh_reader(input_path, mesh_file):
    """
    Read an .STL or .PLY file and return a mesh object
    """

    print("Testing IO for mesh ...")
    mesh = o3d.io.read_triangle_mesh(input_path+mesh_file)
    print(mesh)
    return mesh

def point_cloud_writer(mesh, sample_points):
    """"
    Use Poisson disk sampling to create a point cloud from the mesh
    object using sample points, then write to wdir
    """

    print("Sampling mesh for point cloud ...")
    pcd = mesh.sample_points_poisson_disk(number_of_points=sample_points, init_factor=3)
    # pcd = mesh.sample_points_uniformly(number_of_points=sample_points)
    print(pcd)

    # Paint the point cloud
    pcd.paint_uniform_color([1, 1, 1])
    o3d.io.write_point_cloud("point_cloud.xyz", pcd)
    return pcd
