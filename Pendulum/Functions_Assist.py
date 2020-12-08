#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 14:01:01 2020

@author: vvalera
"""

import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
import sys

# External Functions import
sys.path.append('../../External_Functions/')
from ExternalFunctions import Chi2Regression

# Read in data
def read_data(filename):
    dat = np.genfromtxt(filename, delimiter='\t', names=('n', 't_s'))
    return dat

# Define linear function
def linear_fit(x, p0, p1):
    return p0 + p1*x

# Define a gaussian
def gaussian(x, N, mu, sigma):
    return N * 1.0 / (sigma*np.sqrt(2*np.pi)) * np.exp(-0.5* (x-mu)**2/sigma**2)



# Routine for subplots inside plots
def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height])#,axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax

def FitAndResiduals(filename):
    data_example = read_data(filename)
    n, t = data_example['n'], data_example['t_s']
    if len(n) > 15:
        n = n[:15]
        t = t[:15]
    
    # Perform the fit to a linear function
    chi2_object_lin = Chi2Regression(linear_fit, n, t)
    minuit_lin = Minuit(chi2_object_lin, pedantic=False, p0=5.0, p1=3.0)
    minuit_lin.migrad()
    
    p0 = minuit_lin.values['p0']
    p1 = minuit_lin.values['p1']
    
    sigma_p0 = minuit_lin.errors['p0']
    sigma_p1 = minuit_lin.errors['p1']    
    
    # Compute the residuals
    t_predict = linear_fit(n, *minuit_lin.args)
    residuals = t - t_predict    
    
    return p0, p1, sigma_p0, sigma_p1, residuals

def Fit2Gaussian(bin_centers, counts):
    s_counts = np.sqrt(counts)
    Chi2_object = Chi2Regression(gaussian, bin_centers[counts>0], counts[counts>0], s_counts[counts>0])
    minuit = Minuit(Chi2_object, pedantic=False, N=10, mu=0.0, sigma=0.3, print_level=0) #   
    minuit.migrad()
    N = minuit.values['N']
    mu = minuit.values['mu']
    sigma = minuit.values['sigma']
    return N, mu, sigma, s_counts

def WeightedMean(values, errors):
    if len(values) != len(errors):
        return "Error"
    else:
        a = np.sum(values/errors**2)
        b = np.sum(1.0/errors**2)
        mu = a/b
        sigma = np.sqrt(1.0/b)
    return mu, sigma

