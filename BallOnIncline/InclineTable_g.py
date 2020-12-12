#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 15:52:34 2020

@author: vvalera
"""
import numpy as np
from Acceleration import ComputeFact1
from Theta import GetThetaGonio, GetThetaTrigon

size = 'small'
#size = 'big' #Choose the size of the rolling ball, notice they are two different experiments

Fact1, sigma_Fact1, [sigma2_a, sigma2_D, sigma2_d] = ComputeSi


g = a/(np.sin(theta))*(1 + 2.0/5.0*D**2/(D**2-d**2))
print(g)