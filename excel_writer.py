#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 17:26:31 2020

@author:

    Jack Shine

    Graduate Student Researcher

    National Aerothermochemistry and Hypersonic Labratory

    Texas A & M University
"""
import pandas as pd

def exporter(A, phi, lambda_k, N, M, pts_perturbed, output_path):
    """
    This module exports parameters pertinent to the quasi-random
    rough surface generator
    """

    # Create a Pandas dataframe from the data.
    df_A = pd.DataFrame(A)
    df_phi = pd.DataFrame(phi)
    df_lam = pd.DataFrame({"Wavelength": [lambda_k]})
    df_N = pd.DataFrame({"N": [N]})
    df_M = pd.DataFrame({"M": [M]})
    df_pts = pd.DataFrame(pts_perturbed)

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(output_path+"fourier_parameters.xlsx", engine='xlsxwriter')

    # Convert the dataframe to an XlsxWriter Excel object.
    df_A.to_excel(writer, sheet_name='Amplitude', index=False, header=0)
    df_phi.to_excel(writer, sheet_name='Phase Shift', index=False, header=0)
    df_lam.to_excel(writer, sheet_name='Wavelength', index=False)
    df_N.to_excel(writer, sheet_name='N', index=False)
    df_M.to_excel(writer, sheet_name='M', index=False)
    df_pts.to_excel(writer, sheet_name='Perturbed Points', index=False,header=0)


    # Close the Pandas Excel writer and output the Excel file.
    writer.save()
