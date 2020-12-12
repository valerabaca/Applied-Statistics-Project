#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 12 13:39:10 2020

@author: vvalera
"""

import numpy as np

def WeightedMean(values, errors):
    weights = 1.0/errors**2
    mu = np.average(values, weights=weights)
    sigma = np.sqrt(1.0/np.sum(weights))
    return mu, sigma