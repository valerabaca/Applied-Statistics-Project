#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  7 11:46:51 2020

@author: vvalera
"""
# Standard Imports
import numpy as np
import matplotlib.pyplot as plt
import sys

# External Functions import - Troels Functions
sys.path.append('../External_Functions')
from ExternalFunctions import nice_string_output, add_text_to_ax   # Useful functions to print fit results on figure

# Local external function a.k.a my functions written in Functions_Assist.py (must be in the same folder or add to the sys.path as above)
from Functions_Assist import read_data, linear_fit, add_subplot_axes, FitAndResiduals, Fit2Gaussian, gaussian


# Getting the name of the data files
#Filename1 = '../StopwatchTimer/Measurements_Victor_1.dat' 
Filename1 = '../StopwatchTimer/timer_output.dat' 
Filename2 = '../StopwatchTimer/Measurements_Victor_2.dat'
Filename3 = '../StopwatchTimer/Measurements_Victor_3.dat'

# Here calls the function that makes the fit and compute de residuals, and stores it in a Python Library
DATA1 = {}
DATA1['p0'], DATA1['p1'], DATA1['sigma_p0'], DATA1['sigma_p1'], DATA1['residuals'] = FitAndResiduals(filename=Filename1)

DATA2 = {}
DATA2['p0'], DATA2['p1'], DATA2['sigma_p0'], DATA2['sigma_p1'], DATA2['residuals'] = FitAndResiduals(filename=Filename2)

DATA3 = {}
DATA3['p0'], DATA3['p1'], DATA3['sigma_p0'], DATA3['sigma_p1'], DATA3['residuals'] = FitAndResiduals(filename=Filename3)


# This is dummy step, because the same info could be extracted from the function above, but... whatever
DATA1_read = read_data(Filename1)
n1, t1 = DATA1_read['n'], DATA1_read['t_s']
n1, t1 = n1[:15], t1[:15]

DATA2_read = read_data(Filename2)
n2, t2 = DATA2_read['n'], DATA2_read['t_s']
n2, t2 = n2[:15], t2[:15]

DATA3_read = read_data(Filename3)
n3, t3 = DATA3_read['n'], DATA3_read['t_s']
n3, t3 = n3[:15], t3[:15]
        
# Use the fitted parameters to draw the fitted line usinf the function linear_fit in Functions_Assist.py
n_fit_1 = np.linspace(0.9*n1.min(), 1.1*n1.max())
t_fit_1 = linear_fit(n_fit_1, DATA1['p0'], DATA1['p1'])

n_fit_2 = np.linspace(0.9*n2.min(), 1.1*n2.max())
t_fit_2 = linear_fit(n_fit_2, DATA2['p0'], DATA2['p1'])

n_fit_3 = np.linspace(0.9*n3.min(), 1.1*n3.max())
t_fit_3 = linear_fit(n_fit_3, DATA3['p0'], DATA3['p1'])

# Combine the residuals indo a single Numpy array to make a histogram with them
res_combined = np.concatenate((DATA1['residuals'], DATA2['residuals'], DATA3['residuals']))

# Plotting
# First define a subplot array for the main figures
sig_t = 0.1   # Set your own values...
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(14, 10),
                       gridspec_kw={'height_ratios':[4,1]}, sharex=True)
ax[0].errorbar(n1, t1, yerr=sig_t, color='k', fmt='o')
ax[0].errorbar(n2, t2, yerr=sig_t, color='k', fmt='o')
ax[0].errorbar(n3, t3, yerr=sig_t, color='k', fmt='o')
ax[0].set_xlabel('Timing measurement number')
ax[0].set_ylabel('Time elapsed (s)')
ax[0].set(xlim=(0, 16), ylim=(0, 55))

ax[0].plot(n_fit_1, t_fit_1, '-r')
ax[0].plot(n_fit_2, t_fit_2, '-b')
ax[0].plot(n_fit_3, t_fit_3, '-g')

# Plotting Residuals
ax[1].errorbar(n1, DATA1['residuals'], yerr=sig_t, color='r', fmt='x')
ax[1].axhline(0, xmin=0, xmax=20, color='k')

ax[1].errorbar(n2, DATA2['residuals'], yerr=sig_t, color='b', fmt='x')
ax[1].axhline(0, xmin=0, xmax=20, color='k')

ax[1].errorbar(n3, DATA3['residuals'], yerr=sig_t, color='g', fmt='x')
ax[1].axhline(0, xmin=0, xmax=20, color='k', ls='-')
ax[1].axhline(0+sig_t, xmin=0, xmax=20, color='k', ls='--')
ax[1].axhline(0-sig_t, xmin=0, xmax=20, color='k', ls='--')

# Then we call the function add_subplot_axes in Functions_Assists.py to draw the in-picture figure
subpos = [.65,0.05,0.3,0.4]
subax1 = add_subplot_axes(ax[0],subpos)
subax1.set(title='Residuals Distribution')
counts, bin_edges = np.histogram(res_combined,bins=7)
bin_centers = (bin_edges[1:] + bin_edges[:-1])/2

# Fit the histogram of residuals to a gaussian
N, mu, sigma, s_counts =  Fit2Gaussian(bin_centers=bin_centers, counts=counts)
xaxis = np.linspace(-0.3, 0.3, 1000)
yaxis = gaussian(xaxis, N, mu, sigma)
subax1.plot(xaxis, yaxis, ls='--', color='k')
subax1.errorbar(bin_centers, counts, yerr=s_counts, fmt='x', color='m')

# Adding text to file
d = {r'$T_1$': [DATA1['p1'], DATA1['sigma_p1']],
     r'$T_2$': [DATA2['p1'], DATA2['sigma_p1']],
     r'$T_3$': [DATA3['p1'], DATA3['sigma_p1']],
     }
text = nice_string_output(d, extra_spacing=2, decimals=3)
add_text_to_ax(0.02, 0.97, text, ax[0], fontsize=15)
