#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filename: lc_stack.py
Author: Nat Heddaeus
Date: 2024-11-29
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

# storing variable data as dictionaries of {index, <var>}; it is assumed that indices of desired jd for measurement will be known at this stage
jd = {}
forceddiffimflux = {}
forceddiffimfluxunc = {}
zpdiff = {}
filter = {}

forceddiffimflux_new = {}
flux_by_filter = {}  # {filter, list of flux}
forceddiffimfluxunc_new = {}
unc_by_filter = {}  # {filter, list of unc}


def fill_vars():
    "Method for taking a string of ascii data and converting it to variables."
    for line in data:
        line = line.strip()
        columns = line.split()
        jd[int(columns[0])] = (float(columns[22]))  # creates dictionary of {index, julian day}
        if (columns[24] != 'null'): 
            forceddiffimflux[int(columns[0])] = (float(columns[24]))  # creates dictionary of {index, flux}
        else:
            forceddiffimflux[int(columns[0])] = (None)  # catches unusable data points
        if (columns[25] != 'null'): 
            forceddiffimfluxunc[int(columns[0])] = (float(columns[25]))  # creates dictionary of {index, flux uncertainty}
        else:
            forceddiffimflux[int(columns[0])] = (None)  # catches unusable data points
        zpdiff[int(columns[0])] = (float(columns[20]))  # creates dictionary of {index, zero points}
        filter[int(columns[0])] = (columns[4])  # creates dictionary of {index, filter}


def get_zpavg(start, end):  # start, end inclusive, needs to be updated to calculate by filter
    """Helper method for rescale()."""
    tot = 0
    size = 0
    for index in zpdiff:
        if index in range(start, end + 1):
            if forceddiffimflux[index] != None and forceddiffimfluxunc[index] != None:  
                tot += zpdiff[index]
                size += 1
    return (tot/size)


def rescale(start, end):  # start, end inclusive; need to add index out of bounds errors for < first index in list, > last index in list
    """Given lists of ascii data, make new columns for forcediffimfluxi and forceddiffimfluxunci with rescaled input fluxes and uncertainties."""
    zpavg = get_zpavg(start, end + 1)
    g_list = []  # list of rescaled fluxes in the g band
    r_list = []  # list of rescaled fluxes in the r band
    i_list = []  # list of rescaled fluxes in the i band
    g_unc_list = []  # list of rescaled uncertainties in the g band
    r_unc_list = []  # list of rescaled uncertainties in the r band
    i_unc_list = []  # list of rescaled uncertainties in the i band

    for index in forceddiffimflux:  # check calculations
        if index in range (start, end + 1):
            if forceddiffimflux[index] != None and forceddiffimfluxunc[index] != None:  # skips unusable data points
                new = forceddiffimflux[index]*10**(0.4*(zpavg-zpdiff[index]))
                forceddiffimflux_new[index] = new  # place fluxes on the same photometric zeropoint, may be removable
                # sorts new flux by filter
                if filter[index] == 'ZTF_g':
                    g_list.append(new)
                elif filter[index] == 'ZTF_r':
                    r_list.append(new)
                elif filter[index] == 'ZTF_i':
                    i_list.append(new)
    
    flux_by_filter['ZTF_g'] = g_list
    flux_by_filter['ZTF_r'] = r_list
    flux_by_filter['ZTF_i'] = i_list

    for index in forceddiffimfluxunc:  # check calculations
        if index in range(start, end + 1):
            if forceddiffimfluxunc[index] != None and forceddiffimflux[index] != None:  # skips unusable data points
                new = forceddiffimfluxunc[index]*10**(0.4*(zpavg-zpdiff[index]))
                forceddiffimfluxunc_new[index] = new  # place uncertainties on the same photometric zeropoint, may be removable
                # sorts new uncertainty by filter
                if filter[index] == 'ZTF_g':
                    g_unc_list.append(new)
                elif filter[index] == 'ZTF_r':
                    r_unc_list.append(new)
                elif filter[index] == 'ZTF_i':
                    i_unc_list.append(new)
    
    unc_by_filter['ZTF_g'] = g_unc_list
    unc_by_filter['ZTF_r'] = r_unc_list
    unc_by_filter['ZTF_i'] = i_unc_list


# def collapse_flux():  # will need to be reworked for specific time window; check calculations
#     """Assuming that underlying source is stationary within time window, collapse rescaled single-epoch fluxes using an inverse-variance weighted average."""
#     flux = 0
#     w_tot = 0
#     for index in forceddiffimflux_new:
#         w = 1/(forceddiffimfluxunc_new[index])**2  # weight (uncertainty)
#         w_tot += w  # summation of weighted uncertainties
#         flux += w*forceddiffimflux_new[index]
#     flux = flux / w_tot
#     flux_unc = w_tot**(-1/2)  # effectively forceddiffimfluxunc_new/sqrt(n)
#     return flux, flux_unc


def collapse_flux_by_filter():  # needs testing
    """Assuming that underlying source is stationary within time window, collapse rescaled single-epoch fluxes using an inverse-variance weighted average. Bin Separately by Filter"""
    for filter in flux_by_filter:
        flux = 0
        w_tot = 0
        for point in flux_by_filter[filter]:
            w = 1/(point)**2  # weight (uncertainty)
            w_tot += w
            flux += w*point
        if w_tot != 0:
            flux = flux / w_tot
            flux_unc = w_tot**(-1/2) 
        print(filter + ' flux: ' + str(flux) + ' flux_unc: ' + str(flux_unc))  # will need to convert these to variables

    
# Step 3. Alternative Method for Collapsing Measurements: TBD

# Step 4. Convert Flux to calibrated magnitudes with upper-limits; see section 12 and 13.

def main():
    fill_vars()
    #print(get_zpavg(0, 2335))
    rescale(0, 2335)
    #print(unc_by_filter)
    #print(forceddiffimfluxunc)
    #print(forceddiffimfluxunc_new)
    #print(zpdiff)
    collapse_flux_by_filter()

main()