#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filename: lc_stack.py
Author: Nat Heddaeus
Date: 2024-12-22
Version: 1.0
Description: Given a dataframe with a maxlike light curve, stack the flux. Based on guidelines from "Generating Lightcurves from Forced PSF-fit Photometry on ZTF Difference Images" by Masci et. al, 2022.

Contact: nathedd@unc.edu
"""

import math

# use of user inputs to select desired binning window
file_name = str(input("Enter filename: "))
start = 57 + int(input("Enter starting index (Index will be found in column one of your ascii file, in the rows that do not start with '#'): "))  # starting line 57 appears to be consistent across ZTF files, but may need to be changed in future updates
end = 58 + int(input("Enter ending index: "))  # readlines method appears to be end exclusive, hence 57 + 1
baseline = float(input("If there is any residual baseline (nonzero), input it here. Else, input 0: "))  # you will need to have examined a plot of forceddiffimflux to jd to determine if there is any residual offset in the baseline (Masci et. al section 10)

with open(file_name) as f: # opens .txt file
    data = f.readlines()[start: end]

# storing variable data as dictionaries: key = index, value = column of file_name
jd = {}  # julian day; not used in these methods but may be passed to other files
forceddiffimflux = {}
forceddiffimfluxunc = {}
zpdiff = {}  # zeropoint
filter = {}

forceddiffimflux_new = {}  # not used in these methods but may be passed to other files
flux_by_filter = {}  # {filter, list of flux}
forceddiffimfluxunc_new = {}  # not used in these methods but may be passed to other files
unc_by_filter = {}  # {filter, list of unc}


def fill_vars():  # null value replacement with None may be unneccessary now
    "Takes a string of ascii data and converts it to dictionary variables."
    for line in data:
        line = line.strip()
        columns = line.split()

        jd[int(columns[0])] = (float(columns[22]))  # creates dictionary of {index, julian day}
        
        if (columns[24] != 'null'): 
            forceddiffimflux[int(columns[0])] = (float(columns[24]))  # creates dictionary of {index, flux}
        else:
            forceddiffimflux[int(columns[0])] = (None)  # catches unusable data points; fills with None for consistency between variables
        
        if (columns[25] != 'null'): 
            forceddiffimfluxunc[int(columns[0])] = (float(columns[25]))  # creates dictionary of {index, flux uncertainty}
        else:
            forceddiffimflux[int(columns[0])] = (None)  # catches unusable data points; fills with None for consistency between variables
        
        zpdiff[int(columns[0])] = (float(columns[20]))  # creates dictionary of {index, zero points}
        filter[int(columns[0])] = (columns[4])  # creates dictionary of {index, filter}


def correct_baseline():
    for key in forceddiffimflux:
        """Following an estimate of the baseline level, subtract estimate from differential flux measurements."""
        forceddiffimflux[key] = forceddiffimflux[key] - baseline


def get_zpavg():  # start, end inclusive, needs to be updated to calculate by filter
    """Helper method for rescale()."""
    tot = 0
    size = 0
    for index in zpdiff:
        if forceddiffimflux[index] != None and forceddiffimfluxunc[index] != None:  
            tot += zpdiff[index]  # zpdiff appears to be more or less consistent across filters
            size += 1
    return (tot/size)


def rescale():  # start, end inclusive; need to add index out of bounds errors for < first index in list, > last index in list
    """Given lists of ascii data, make new columns for forcediffimfluxi and forceddiffimfluxunci with rescaled input fluxes and uncertainties."""
    zpavg =  26.13802754993896
    g_list = []  # list of rescaled fluxes in the g band
    r_list = []  # list of rescaled fluxes in the r band
    i_list = []  # list of rescaled fluxes in the i band
    g_unc_list = []  # list of rescaled uncertainties in the g band
    r_unc_list = []  # list of rescaled uncertainties in the r band
    i_unc_list = []  # list of rescaled uncertainties in the i band

    for index in forceddiffimflux:  # check calculations
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


def collapse_flux_by_filter():  # needs testing
    """Assuming that underlying source is stationary within time window, collapse rescaled single-epoch fluxes using an inverse-variance weighted average. Bin Separately by Filter"""
    for filter in flux_by_filter:
        flux = 0
        w_tot = 0
        idx = 0
        
        for point in flux_by_filter[filter]:
            w = 1/((unc_by_filter[filter][idx])**2)
            flux += w*point
            w_tot += w
            idx += 1

        if w_tot != 0:
            flux = flux / w_tot
            flux_unc = w_tot**(-1/2) 
            print(filter + ' flux: ' + str(flux) + ' flux_unc: ' + str(flux_unc))  # will need to convert these to variables


def cal_mag():  # will need to use the zpavg value assumed when rescaling the input fluxes
    """Obtaining calibrated magnitudes (for transients)."""
    if (('forceddiffimflux' / 'forceddiffimfluxunc') > 5):  # will need to take flux and unc from collapse_flux_by_filter; 5 is the signal to noise threshold for declaring a measurement a "non-detection", so that it can be assigned an upper-limit (see Masci et. al)
        # confident detection, plot magnitude with error bars
        mag = 'zpdiff' - 2.5*log10['forceddiffimflux']
        sigma = 1.0857 * 'forceddiffimfluxunc'/'forceddiffimflux'
    else:
        # compute upper flux limits and plot as arrow
        mag = zp = 2.5*log10[3 * 'forceddiffimfluxunc']  # 3 is the actual signal to noise ratio to use when computing SNU-sigma upper-limit


def main():
    fill_vars()
    if baseline != 0:
        correct_baseline()
    rescale()
    collapse_flux_by_filter()

main()