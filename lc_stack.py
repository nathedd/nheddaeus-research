#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filename: lc_stack.py
Author: Nat Heddaeus
Date: 2024-10-31
Version: 1.0
Description: Given a dataframe with a maxlike light curve, stack the flux. Based on guidelines from "Generating Lightcurves from Forced PSF-fit Photometry on ZTF Difference Images" by Masci et. al, 2022.

Contact: nathedd@unc.edu
"""


from astropy.io import ascii  # ZTF data file will be in ascii format
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt


with open('results_ZTF24aapvieu.txt') as f: # open .txt file
    data = f.readlines()[57:2393]  # readlines may be a temporary measure in place of a better method of reading the file; index 57 may be specific to this file

# storing variable data as dictionaries of {index, <var>} for now
jd = {}
forceddiffimflux = {}
forceddiffimfluxunc = {}
zpdiff = {}
filter = {}


def fill_vars(string):
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

fill_vars(data)
# Step 1. Place fluxes on same standard zpavg (a value close to zpdiffi values); rescale input fluxes and uncertainties
# forceddiffimflux_newi = forceddiffimfluxi*10**(0.4(zpavg-zpdiffi))
# forceddiffimfluxunc_newi = forceddiffimfluxunci*10**(0.4(zpavg-zpdiffi))
#def rescale(data):
   # """Given an astropy ascii table, replace columns for forcediffimfluxi and forceddiffimfluxunci with rescaled input fluxes and uncertainties."""

# Step 2. Combine Measurements: Assume underlying source is stationary w/in a time window; collapse fluxes (step 1) using an inverse-variance weighted average
# Flux = summation (wiforceddiffimflux_newi) / summation (wi)
# wi = 1/(diffimfluxunc_newi)**2
# Fluxunc = [summation(wi)]^-1/2

# Step 3. Alternative Method for Collapsing Measurements: TBD

# Step 4. Convert Flux to calibrated magnitudes with upper-limits; see section 12 and 13.