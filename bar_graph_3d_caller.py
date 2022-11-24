#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 20:33:54 2022

@author: Jack
"""

# import tkinter as tk
# from tkinter import filedialog

# root = tk.Tk()
# root.withdraw()

# file_path = filedialog.askopenfilename()

import pandas as pd
from bar_graph_3d import bar_graph

file_path = '/Volumes/GoogleDrive/My Drive/Python/large files/2.0x Perturbed Mesh/fourier_parameters.xlsx'

A = pd.read_excel(file_path, sheet_name='Amplitude',header=None).to_numpy()
phi = pd.read_excel(file_path, sheet_name='Phase Shift',header=None).to_numpy()
lambda_k = pd.read_excel(file_path, sheet_name='Wavelength',header=0).astype(float).iloc[0,0]
N = pd.read_excel(file_path, sheet_name='N',header=0).astype(float).iloc[0,0]
M = pd.read_excel(file_path, sheet_name='M',header=0).astype(float).iloc[0,0]

# Plot the A matrix
bar_graph(A, N, M, lambda_k)

