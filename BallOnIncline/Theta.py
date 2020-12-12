#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 15:49:04 2020

@author: vvalera
"""

import numpy as np
from FunctionsAssist import WeightedMean

Trigon = {}
Trigon['Victor'] = [88.89,22.20] # Base, Height cm
Trigon['Peiyuan'] = [88.90,22.20] # cm
Trigon['Sippo'] = [88.90,22.30] # cm
err = 0.05 #cm

Gonio_small_table={}
Gonio_small_table['Victor'] = [88.00, 91.00]
Gonio_small_table['Peiyuan'] = [89.00, 90.00]
Gonio_small_table['Sippo'] = [88.50, 91.00]

Gonio_small_slope = {}
Gonio_small_slope['Victor'] = [78.0, 105.0]
Gonio_small_slope['Peiyuan'] = [78.0, 105.0]
Gonio_small_slope['Sippo'] = [78.0, 105.0]

Gonio_large_table={}
Gonio_large_table['Victor'] = [89.80, 90.10]
Gonio_large_table['Peiyuan'] = [90.00, 90.20]
Gonio_large_table['Sippo'] = [89.50, 90.00]

Gonio_large_slope = {}
Gonio_large_slope['Victor'] = [75.9, 103.9]
Gonio_large_slope['Peiyuan'] = [75.9, 103.9]
Gonio_large_slope['Sippo'] = [75.9, 104.0]


Err_small = 0.5
Err_large = 0.05


def GetThetaTrigon(name):
    theta = np.arctan(Trigon[name][1]/Trigon[name][0])
    sigma_theta = 1.0/np.sqrt(Trigon[name][0]**2 + Trigon[name][1]**2) * err
    return theta, sigma_theta

def GetThetaGonio(name):
    theta_small = (Gonio_small_slope[name][0] + (180.0 - Gonio_small_slope[name][1]))*0.5
    sigma_theta_small = np.sqrt(2)*Err_small
    theta_large = (Gonio_large_slope[name][0] + (180.0 - Gonio_large_slope[name][1]))*0.5
    sigma_theta_large = np.sqrt(2)*Err_large
    delta_theta_small = (Gonio_small_table[name][0] + (180.0 - Gonio_small_table[name][1]))*0.5
    delta_theta_large = (Gonio_large_table[name][0] + (180.0 - Gonio_large_table[name][1]))*0.5
    
    ave_theta, sigma_theta = WeightedMean(np.array([theta_small, theta_large]), np.array([sigma_theta_small, sigma_theta_large]))
    ave_delta_theta, sigma_delta_theta =WeightedMean(np.array([delta_theta_small, delta_theta_large]), np.array([sigma_theta_small, sigma_theta_large]))
    theta = 90.0 - ave_theta
    delta_theta = 90.0 - ave_delta_theta
    
    # convert to rad
    
    theta = theta*np.pi/180.0
    delta_theta = delta_theta*np.pi/180.0
    sigma_theta = sigma_theta*np.pi/180.0
    sigma_delta_theta = sigma_delta_theta*np.pi/180.0
    
    return theta, sigma_theta, delta_theta, sigma_delta_theta


