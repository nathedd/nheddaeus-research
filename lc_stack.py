#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filename: lc_stack.py
Author: Nat Heddaeus
Date: 2024-11-05
Version: 1.0
Description: Given a dataframe with a maxlike light curve, stack the flux. Based on guidelines from "Generating Lightcurves from Forced PSF-fit Photometry on ZTF Difference Images" by Masci et. al, 2022.

Contact: nathedd@unc.edu
"""


from astropy.io import ascii  # ZTF data file will be in ascii format
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt


with open('results_ZTF24aapvieu.txt') as f: # open .txt file; need to find means of inputting file
    data = f.readlines()[57:2393]  # readlines may be a temporary measure in place of a better method of reading the file; index 57 may be specific to this file

# storing variable data as dictionaries of {index, <var>} for now; may want to change indices to jd
jd = {}
forceddiffimflux = {}
forceddiffimfluxunc = {}
zpdiff = {}
filter = {}

forceddiffimflux_new = {}
forceddiffimfluxunc_new = {}


def fill_vars():
    "Method for taking a string of ascii data and converting it to variables."
    for line in data:
        line = line.strip()
        columns = line.split()
        jd[int(columns[0])] = (float(columns[22]))
        if (columns[24] != 'null'):  # this may be a temporary measure, assuming null values will not be included at this stage
            forceddiffimflux[int(columns[0])] = (float(columns[24]))
        else:
            forceddiffimflux[int(columns[0])] = (None)
        if (columns[25] != 'null'):
            forceddiffimfluxunc[int(columns[0])] = (float(columns[25]))  # this may be a temporary measure, assumming null values will not be included at this stage
        else:
            forceddiffimflux[int(columns[0])] = (None)
        zpdiff[int(columns[0])] = (float(columns[20]))
        filter[int(columns[0])] = (columns[4])


def get_zpavg():
    """Helper method for rescale()."""
    tot = 0
    size = 0
    for index in zpdiff:
        if forceddiffimflux[index] != None:  # if statement will probably be unnecessary later on
            tot += zpdiff[index]
            size += 1
    return (tot/size)


def rescale():  # need to add index out of bounds errors for < first index in list, > last index in list
    """Given lists of ascii data, make new columns for forcediffimfluxi and forceddiffimfluxunci with rescaled input fluxes and uncertainties."""
    zpavg = get_zpavg()
    for index in forceddiffimflux:  # check calculations
        if forceddiffimflux[index] != None:  # if statement will probably be unnecessary later on
            forceddiffimflux_new[index] = forceddiffimflux[index]*10**(0.4*(zpavg-zpdiff[index]))  # place fluxes on the same photometric zeropoint
    
    for index in forceddiffimfluxunc:  # check calculations
        if forceddiffimfluxunc[index] != None:  # if statement will probably be unnecessary later on
            forceddiffimfluxunc_new[index] = forceddiffimfluxunc[index]*10**(0.4*(zpavg-zpdiff[index]))  # place uncertainties on the same photometric zeropoint


def collapse_flux():  # will need to be reworked for specific time window; check calculations
    """Assuming that underlying source is stationary within time window, collapse rescaled single-epoch fluxes using an inverse-variance weighted average."""
    flux = 0
    w_tot = 0
    for index in forceddiffimflux_new:
        w = 1/(forceddiffimfluxunc_new[index])**2  # weight (uncertainty)
        w_tot += w  # summation of weighted uncertainties
        flux += w*forceddiffimflux_new[index]
    flux = flux / w_tot
    flux_unc = w_tot**(-1/2)  # effectively forceddiffimfluxunc_new/sqrt(n)
    return flux, flux_unc

# Step 3. Alternative Method for Collapsing Measurements: TBD

# Step 4. Convert Flux to calibrated magnitudes with upper-limits; see section 12 and 13.

def main():
    fill_vars(data)
    #print(get_zpavg())
    rescale()
    #print(forceddiffimfluxunc)
    #print(forceddiffimfluxunc_new)
    #print(zpdiff)
    collapse_flux()

main()