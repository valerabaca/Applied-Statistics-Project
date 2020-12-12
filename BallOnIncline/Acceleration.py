#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 10:20:16 2020

@author: vvalera
"""
import numpy as np
from scipy import odr
from FunctionsAssist import WeightedMean

Victor = [19.71, 32.78, 45.8, 58.85, 71.85]  # cm
Peiyuan = [19.75, 32.79, 45.82, 58.83, 71.85] # cm
Sippo = [19.75, 32.85, 45.75, 58.85, 71.75] # cm
Positions = {'Victor':Victor, 'Peiyuan':Peiyuan, 'Sippo': Sippo}
x_err = 0.05*np.ones_like(Victor) # cm

D_small, sigma_D_small = 11.9/10.0, 0.1/10.0 # cm
D_big, sigma_D_big = 18.8/10.0, 0.1/10.0 # cm

d = {'Victor': [6.0, 5.6],'Peiyuan': [6.0,6.0], 'Sippo': [5.9, 6.0]}
sigma_d = np.ones(2)*0.1/10.0 #cm

def quad_fit(p,t):
    p0, p1, p2 = p
    return p0 + p1*t + 0.5*p2*t**2

def read_csv(filename):
    """Read CSV from Waveforms"""
    dat = np.genfromtxt(filename, delimiter=',', skip_header=13, names=True)
    time = dat['Time_s']
    voltage = dat['Channel_1_V']
    return time, voltage

def find_midpoints(filename, show_plot=True):
    """Find midpoints -- quick and dirty"""
    time, voltage = read_csv(filename)
    T = [[],[],[],[],[]]
    i = 0
    eps = 0
    for v, t in zip(voltage, time):
        if v < 4.7 and eps == 1:
            eps = 0
            i+=1
        if v > 4.7:
            T[i].append(t)
            if eps == 0:
                eps = 1
    t_peak = np.zeros(5)
    t_peak_sigma = np.zeros(5)
    i = 0
    for peak in T:
        t_peak[i] = np.mean(peak)
        t_peak_sigma[i] = (peak[-1] - peak[0])/2
        i+=1
    return t_peak, t_peak_sigma


def FitAndResiduals(t,t_sigma, x, x_sigma):  
    quad_model = odr.Model(quad_fit)
    data = odr.RealData(t, x, sx=t_sigma, sy=x_sigma)
    Odr = odr.ODR(data, quad_model, beta0=[0., 1., 1.])
    out = Odr.run()
    popt = out.beta
    perr = out.sd_beta
    a = popt[2]
    sigma_a = perr[2]
    return a, sigma_a

def ComputeA(name,size):
    A = np.zeros(3)
    sigma_A = np.zeros_like(A)
    for i in range(3):
        filename = 'DATA/measure' + str(i+1) + '_' + size + '.csv'
        time, voltage = read_csv(filename)
        timepass, timepass_sig = find_midpoints(filename)
        a, sigma_a = FitAndResiduals(timepass, timepass_sig, Positions[name], x_err) # cm/s^2
        A[i] = a
        sigma_A[i] = sigma_a
    a, sigma_a = WeightedMean(A, sigma_A)
    return a, sigma_a

def ComputeFact1(name, size):
    # a*[1 + 2/5*D**2/(D**2+d**2)]
    a, sigma_a = ComputeA(name, size)
    d_name = np.array(d[name])
    d_name,sigma_d_name = WeightedMean(d_name, sigma_d)
    if size == 'small':
        D, sigma_D = D_small, sigma_D_small
    if size == 'big':
        D, sigma_D = D_big, sigma_D_big
    Fact1 = a*(1.0 + 2.0/5.0*D**2/(D**2+d_name**2))   
    sigma2_a = (1.0 + 2.0*D**2/(5.0*(D**2-d_name**2))) * sigma_a**2 # A term 1/sin(\theta + \Delta \theta)**2 is missing
    sigma2_D = a**2.0*((4*D/(5*(D**2-d_name**2))) - (4*D**3/(5*(D**2-d_name**2))))**2*sigma_D**2 # A term 1/sin(\theta + \Delta \theta)**2 is missing
    sigma2_d = 16*a**2*D**4*d_name**2/(25*(D**2 - d_name**2)**4)*sigma_d_name**2
    sigma2_Fact1 = sigma2_a + sigma2_D  + sigma2_d
    sigma_Fact1 = np.sqrt(sigma2_Fact1)
    return Fact1, sigma_Fact1, [sigma2_a, sigma2_D, sigma2_d]


name = 'Victor'
print(ComputeFact1(name, 'small'))
print(ComputeFact1(name, 'big'))