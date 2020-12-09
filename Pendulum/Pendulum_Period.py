#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 11:46:51 2020

@author: vvalera
"""
# Standard Imports
import numpy as np
# Local external function a.k.a my functions written in Functions_Assist.py (must be in the same folder or add to the sys.path as above)
from Functions_Assist import FitAndResiduals, WeightedMean, Chi2Test

def ComputeT(name):    
    # Getting the name of the data files
    Filename1 = '../StopwatchTimer/Measurements_' + name + '_1.dat' 
    Filename2 = '../StopwatchTimer/Measurements_' + name + '_2.dat' 
    Filename3 = '../StopwatchTimer/Measurements_' + name + '_3.dat'
    
    DATA1, DATA2, DATA3 = {},{},{}
    
    # Here calls the function that makes the fit and compute de residuals, and stores it in a Python Library
    DATA1['p0'], DATA1['p1'], DATA1['sigma_p0'], DATA1['sigma_p1'], DATA1['residuals'], DATA1['Chi2'], DATA1['Prob'] = FitAndResiduals(filename=Filename1)    
    DATA2['p0'], DATA2['p1'], DATA2['sigma_p0'], DATA2['sigma_p1'], DATA2['residuals'], DATA2['Chi2'], DATA2['Prob'] = FitAndResiduals(filename=Filename2)
    DATA3['p0'], DATA3['p1'], DATA3['sigma_p0'], DATA3['sigma_p1'], DATA3['residuals'], DATA3['Chi2'], DATA3['Prob'] = FitAndResiduals(filename=Filename3)
    
    t = np.array([DATA1['p1'], DATA2['p1'], DATA3['p1']])
    sigma_t = np.array([DATA1['sigma_p1'], DATA2['sigma_p1'], DATA3['sigma_p1']])
    
    T, sigma_T = WeightedMean(t, sigma_t)
    #T_test = Chi2Test(t, sigma_t)
    
    return T, sigma_T