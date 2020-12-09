#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 16:02:15 2020

@author: vvalera
"""
import numpy as np
from iminuit import Minuit
from scipy import stats
import sys
# External Functions import
from Functions_Assist import constant, LowStatsSTD

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

hook = {'Victor': 1.140, 'Sippo': 1.160, 'Peiyuan': 1.165}
weight = {'Victor': 3.460, 'Sippo': 3.450, 'Peiyuan': 3.445}

Hook = np.array(list(hook.values()))
Weight = np.array(list(weight.values()))

h_err = LowStatsSTD(Hook)
w_err = LowStatsSTD(Weight)

def FitConstant(values, errors):
    Npoints = len(values)
    x = np.arange(Npoints)
    p_0 = np.mean(values)
    chi2_object = Chi2Regression(constant, x, values, errors)
    minuit_cte = Minuit(chi2_object, pedantic=False, p0 = p_0)
    minuit_cte.migrad()
    mean = minuit_cte.values['p0']
    sigma = minuit_cte.errors['p0']
    return mean, sigma

def ComputeL(name):
    Filename = 'DataLength/' + str(name) + '_Measurements.dat'
    l1, l2 = read_data(Filename)[:,0], read_data(Filename)[:,1]
    l1_err, l2_err = np.ones_like(l1)*LowStatsSTD(l1), np.ones_like(l2)*LowStatsSTD(l1)
    L1, sigma_L1 = FitConstant(l1, l1_err)
    L2, sigma_L2 = FitConstant(l2, l2_err)
    
    w = weight[name]
    h = hook[name]    
    
    L = L1 + L2 + h + 0.5*w
    
    sigma2_L = sigma_L1**2 + sigma_L2**2 + h_err**2 + (0.5*w_err)**2
    sigma_L = np.sqrt(sigma2_L)
    
    return L, sigma_L


