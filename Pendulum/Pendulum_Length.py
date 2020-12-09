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
from Functions_Assist import WeightedMean, LowStatsSTD
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

DATA1['l'] = np.mean(DATA1['l1']) + np.mean(DATA1['l2'])
DATA1['l_err'] = np.sqrt(2)*l_err

DATA2['l'] = np.mean(DATA2['l1']) + np.mean(DATA2['l2'])
DATA2['l_err'] = np.sqrt(2)*l_err

DATA3['l'] = np.mean(DATA3['l1']) + np.mean(DATA3['l2'])
DATA3['l_err'] = np.sqrt(2)*l_err

DATA1['h'], DATA1['h_err'] = 1.140, 0.005
DATA2['h'], DATA2['h_err'] = 1.160, 0.005
DATA3['h'], DATA3['h_err'] = 1.165, 0.005

DATA1['w'], DATA1['w_err'] = 3.460, 0.005
DATA2['w'], DATA2['w_err'] = 3.450, 0.005
DATA3['w'], DATA3['w_err'] = 3.445, 0.005

# L = L1 + L2 + hook[0] + w[0]/2

l, sigma_l = WeightedMean(np.array([DATA1['l'], DATA2['l'], DATA3['l']]),
                          np.array([DATA1['l_err'], DATA2['l_err'], DATA3['l_err']]))
h, sigma_h = WeightedMean(np.array([DATA1['h'], DATA2['h'], DATA3['h']]),
                          np.array([DATA1['h_err'], DATA2['h_err'], DATA3['h_err']]))
w, sigma_w = WeightedMean(np.array([DATA1['w'], DATA2['w'], DATA3['w']]),
                          np.array([DATA1['w_err'], DATA2['w_err'], DATA3['w_err']]))

L = l + h + 0.5*w
sigma_L = np.sqrt(sigma_l**2 + sigma_h**2 + (0.5*sigma_w)**2)

t = np.array([2.761327857143244, 2.7591366917295193, 2.7684828571429394,
              2.7660266253872243, 2.7625892481204497, 2.7612087719298453,
              2.772606428571263, 2.778663333333263, 2.7570821078432926])
t_err= np.array([0.059761432806335064, 0.038778331473693695, 0.036037498315646795,
                 0.04543108453583409, 0.03877833763092251, 0.041885391948211156,
                 0.059761429798607296, 0.041885390717901254, 0.04950737707058824])
Period = WeightedMean(t, t_err)

T = Period[0]
sigma_T = Period[1]

g = L*(2*np.pi/T)**2/100.0 # From cm to m

sigma2_g = (2*np.pi/T)**4*sigma_L**2 + (2*L*(2*np.pi)**2/T**3)**2*sigma_T**2
sigma_g = np.sqrt(sigma2_g)/100.0
print(g, sigma_g)
