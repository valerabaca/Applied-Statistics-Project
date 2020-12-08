#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 16:02:15 2020

@author: vvalera
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
# External Functions import
from Functions_Assist import WeightedMean
sys.path.append('../External_Functions')
from ExternalFunctions import Chi2Regression

# Read in data
def read_data(filename):
    dat = np.genfromtxt(filename, delimiter='\t')
    return dat

# Defining the systematic erros

l_err = 0.05 # Error with the ruler 0.05 cm
w_err = 0.05/10.0 # Error with the caliper 0.05 mm (converted to cm)
h_err = 0.05/10.0 # Error with the caliper 0.05 mm (converted to cm)

Filename1 = 'DataLength/Victor_Measurements.dat'
Filename2 = 'DataLength/Sippo_Measurements.dat'
Filename3 = 'DataLength/Peiyuan_Measurements.dat'

DATA1, DATA2, DATA3 = {}, {}, {}

DATA1['l1'], DATA1['l2'] = read_data(Filename1)[:,0], read_data(Filename1)[:,1]
DATA2['l1'], DATA2['l2'] = read_data(Filename2)[:,0], read_data(Filename2)[:,1]
DATA3['l1'], DATA3['l2'] = read_data(Filename3)[:,0], read_data(Filename3)[:,1]

L1 = np.average(l1)
L2 = np.average(l2)
l1_err = np.std(l1)
l2_err = np.std(l2)

w = np.array([34.50, 35,35, 34.60])/10.0 #to cm
hook = np.array([11.60, 11.40, 11.65])/10.0 #to cm

w_err = 0.05/10.0 #to cm
hook_err = 0.05/10.0 #to cm 

L = L1 + L2 + hook[0] + w[0]/2

g = L*(2*np.pi/2.761)**2/100.0 # From cm to m
print(g)