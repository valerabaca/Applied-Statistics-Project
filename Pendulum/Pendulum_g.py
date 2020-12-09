#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 16:52:54 2020

@author: vvalera
"""
import numpy as np
from Pendulum_Length import ComputeL
from Pendulum_Period import ComputeT
from Functions_Assist import WeightedMean, Chi2Test

L, T, sigma_L, sigma_T= {}, {}, {}, {}
names = ['Victor', 'Sippo', 'Peiyuan']

for name in names:
    L[name], sigma_L[name] = ComputeL(name)
    T[name], sigma_T[name] = ComputeT(name)

L_value, L_sigma = np.array(list(L.values())), np.array(list(sigma_L.values()))
T_value, T_sigma = np.array(list(T.values())), np.array(list(sigma_T.values()))

L_test = Chi2Test(L_value, L_sigma)
T_test = Chi2Test(T_value, T_sigma)

L, sigma_L = WeightedMean(L_value, L_sigma)
T, sigma_T = WeightedMean(T_value, T_sigma)

g = L*(2*np.pi/T)**2/100.0 # From cm to m

sigma2_g = (2*np.pi/T)**4*sigma_L**2 + (2*L*(2*np.pi)**2/T**3)**2*sigma_T**2
sigma_g = np.sqrt(sigma2_g)/100.0

print("L = {:.3f} +- {:.3f}, \t Chi2 = {:.3f} and p = {:.3f}".format(L, sigma_L, L_test[0], L_test[1]))
print("T = {:.3f} +- {:.3f}, \t Chi2 = {:.3f} and p = {:.3f}".format(T, sigma_T, T_test[0], T_test[1]))
print("g = {:.3f} +- {:.3f}".format(g, sigma_g))